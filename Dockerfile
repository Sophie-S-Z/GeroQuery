# Reproducible container for the GeroQuery API.
FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    GEROQUERY_DATA_HOME=/data

WORKDIR /app

# Install dependencies first for better layer caching.
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Install the package.
COPY pyproject.toml README.md ./
COPY geroquery ./geroquery
RUN pip install --no-cache-dir .

# Materialize the bundled data slice at build time so first request is fast.
RUN python -c "from geroquery.store import GeroStore; GeroStore().build()"

EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
  CMD python -c "import httpx,sys; sys.exit(0 if httpx.get('http://localhost:8000/healthz').status_code==200 else 1)"

CMD ["uvicorn", "geroquery.api.app:app", "--host", "0.0.0.0", "--port", "8000"]
