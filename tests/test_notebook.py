"""The tutorial notebook is the documented end-to-end tour, and nothing ran it.

Its code cells are executed in order, offline, against the bundled store, so an
API change that breaks the tour fails here instead of in a reader's hands. The
test also pins what the tour shows: its aging-clock cell must use a calibrated
clock. It used to apply the hand-set reference demo to real NHANES adults and
print a mean age acceleration of +46 years as though it meant something.
"""

from __future__ import annotations

import json
from pathlib import Path

NOTEBOOK = Path(__file__).resolve().parents[1] / "notebooks" / "01_end_to_end.ipynb"


def test_the_tutorial_notebook_runs_and_its_clock_is_calibrated():
    cells = json.loads(NOTEBOOK.read_text(encoding="utf-8"))["cells"]
    namespace: dict = {}
    for index, cell in enumerate(cells):
        if cell["cell_type"] == "code":
            source = "".join(cell["source"])
            exec(compile(source, f"{NOTEBOOK.name}[{index}]", "exec"), namespace)

    # A calibrated clock averages near zero acceleration on the population it
    # describes; an uncalibrated one is off by decades.
    assert abs(namespace["res"]["mean_age_acceleration"]) < 5
