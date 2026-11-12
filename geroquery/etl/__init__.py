"""Offline, batch ETL that builds the data layer. Never on the request path.

Each step is a ``python -m geroquery.etl.<module>`` target; the Makefile's
``data``, ``signatures``, ``idmap``, ``crosslayer`` and ``frontend-data`` targets
run them in order.
"""
