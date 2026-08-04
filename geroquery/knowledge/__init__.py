"""Citation and vocabulary layer over the ingested data.

This module used to be GeroQuery's *source of truth*: `aging_knowledge.py` held
hand-written, hand-cited per-gene evidence, because the environment it was
authored in had outbound access to biomedical APIs blocked, so real data could
not be fetched.

It is no longer the source of truth. Gene-level evidence now comes from the
checksum-pinned GEO DataSets panel (`sources/geo.py`) and the HAGR releases
(`sources/hagr.py`), which are measured and re-derivable rather than curated by
hand. What survives here is the part the ingestion does *not* provide:

* :data:`REFERENCES` — real, verifiable literature citations, used to point a
  reader at the papers behind a method or a database. Never an invented PMID:
  where the identifier is not certain, the link is a title search that resolves
  to the same paper.
* :data:`HALLMARKS` — the López-Otín hallmarks-of-aging taxonomy. A published
  vocabulary, not data, and the natural way to group findings for a reader.

**Deliberately not carried over:** the per-gene `KNOWLEDGE` table and its
`Evidence` records. Those asserted directions and strengths for individual genes
by hand, which is exactly what the ingested panel now measures — and disagrees
with in at least one prominent case (see `docs/RESULTS_GEO_SIGNATURES.md` on
CDKN2A/p16). Keeping both would mean shipping two answers to the same question,
one of them unfalsifiable.
"""

from __future__ import annotations

from .hallmarks import HALLMARKS
from .references import REFERENCES, Reference, ref

__all__ = ["HALLMARKS", "REFERENCES", "Reference", "ref"]
