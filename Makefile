.PHONY: install dev data api dashboard test lint format typecheck ci clean

install:
	pip install -e .

dev:
	pip install -e ".[dev,ui]"
	pre-commit install || true

data:  ## regenerate the bundled harmonized slice + build the store
	python -m geroquery.etl.build_fixtures
	python -c "from geroquery.store import GeroStore; print('data version', GeroStore().build().version())"

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
