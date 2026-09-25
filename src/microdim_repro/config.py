from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json

PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = PACKAGE_DIR / "datasets.json"

@dataclass(frozen=True)
class Cohort:
    id: str
    label: str
    N: int
    group_column: str


def load_config(path: str | Path | None = None) -> dict:
    path = Path(path) if path is not None else DEFAULT_CONFIG
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def cohorts(path: str | Path | None = None) -> list[Cohort]:
    cfg = load_config(path)
    return [Cohort(**c) for c in cfg["cohorts"]]
