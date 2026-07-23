"""GeroQuery — open-source multi-omic & clinical aging data aggregator.

The package is organized around eight deep modules (Ousterhout's sense): each
hides substantial complexity behind a small, stable interface and is testable in
isolation.

    ui  ->  api  ->  (store, clocks, resilience, harmonize)  ->  (sources, idmap)

Lower modules never import higher ones.
"""

__version__ = "1.0.0"

# Semantic version of the *harmonized data layer*. Bump whenever the ETL that
# regenerates the Parquet/relational store changes in a way that alters results.
DATA_VERSION = "2024.1"

__all__ = ["__version__", "DATA_VERSION"]
