"""M7 api — FastAPI service layer.

Import `geroquery.api.app` for the HTTP app and `geroquery.api.service` for the
service. Nothing is imported here: the dashboard and the MCP tools use only the
service, and an eager `from .app import ...` built a FastAPI app in every one of
them (and crashed the dashboard container when streamlit pulled in a starlette
that FastAPI 0.111 cannot use).
"""
