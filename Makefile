.PHONY: install dev data data-offline crosslayer fetch idmap signatures api dashboard frontend-data frontend test lint format typecheck ci clean

install:
	pip install -e .

dev:
	pip install -e ".[dev,ui]"
	pre-commit install || true

fetch:  ## download + SHA-256 verify every pinned upstream artifact
	GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.fetch_artifacts

idmap:  ## regenerate the bundled gene + AnAge tables from real sources
	GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.build_idmap

signatures:  ## GEO aging panel + HAGR -> signatures, studies, curated knowledge
	GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.build_signatures

crosslayer:  ## NHANES 1999-2002 DNAm clocks + health state + mortality, joined
	GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.build_crosslayer

data: fetch  ## fetch every real upstream, rebuild every table, build the store
	GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.build_data
	GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.build_signatures
	GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.build_idmap
	GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.build_crosslayer
	python -m geroquery.etl.build_fixtures
	python -c "from geroquery.store import GeroStore; print('data version', GeroStore().build().version())"

data-offline:  ## build the store with no network, from the committed real samples
	python -m geroquery.etl.build_fixtures
	python -c "from geroquery.store import GeroStore; print('data version', GeroStore().build().version())"

frontend-data:  ## export the static Parquet/JSON the browser app reads
	python -m geroquery.etl.build_frontend_data

frontend: frontend-data  ## build the static site into frontend/dist
	cd frontend && npm install && npm run build

api:  ## run the FastAPI service (http://localhost:8000/docs)
	uvicorn geroquery.api.app:app --reload --port 8000

dashboard:  ## run the Streamlit dashboard
	streamlit run geroquery/ui/streamlit_app.py

test:
	pytest --cov=geroquery --cov-report=term-missing

lint:
	ruff check geroquery tests

format:
	black geroquery tests
	ruff check --fix geroquery tests

typecheck:
	mypy geroquery

ci: lint typecheck test

clean:
	rm -rf geroquery/_build geroquery/_cache .pytest_cache .ruff_cache .mypy_cache .coverage
