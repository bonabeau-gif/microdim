from __future__ import annotations
from pathlib import Path
from urllib.request import urlretrieve
import pandas as pd
import numpy as np

from .config import load_config, cohorts


def download_file(url: str, destination: str | Path, overwrite: bool = False) -> Path:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        return destination
    urlretrieve(url, destination)
    return destination


def download_public_data(data_dir: str | Path, config_path=None, overwrite=False, include_cohorts=True, include_micom=True):
    """Download only the public files needed by the manuscript analyses.

    The full curated repositories are not cloned. For the eight cohorts this
    retrieves genus abundances, metabolite tables, and metadata. For the MICOM
    proof-of-concept this retrieves the public exchange-flux and diversity files.
    """
    data_dir = Path(data_dir)
    cfg = load_config(config_path)
    downloaded = []
    if include_cohorts:
        repo = cfg["cohort_repository"]
        branch = cfg.get("cohort_branch", "main")
        for c in cohorts(config_path):
            base = f"https://raw.githubusercontent.com/{repo}/{branch}/data/processed_data/{c.id}"
            cdir = data_dir / "cohorts" / c.id
            for name in ("genera.tsv", "mtb.tsv", "metadata.tsv"):
                downloaded.append(download_file(f"{base}/{name}", cdir / name, overwrite))
    if include_micom:
        mdir = data_dir / "micom"
        downloaded.append(download_file(cfg["micom_exchange_url"], mdir / "exchanges_d0_high_fiber.csv", overwrite))
        downloaded.append(download_file(cfg["micom_diversity_url"], mdir / "vendor_diversity.csv", overwrite))
    return downloaded


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    # Curated cohort files are tab-delimited. sep=None is intentionally avoided
    # because large tables make dialect inference slow and fragile.
    return pd.read_csv(path, sep="\t", index_col=0, low_memory=False)


def orient_samples_by_features(table: pd.DataFrame, sample_ids) -> pd.DataFrame:
    """Orient a matrix as samples x features using metadata sample IDs.

    The curated collection contains tables written in several conventions.
    This orientation detector makes the Python reproduction independent of
    whether samples occupy rows or columns in a particular source file.
    """
    ids = set(map(str, sample_ids))
    row_hits = sum(str(x) in ids for x in table.index)
    col_hits = sum(str(x) in ids for x in table.columns)
    if col_hits > row_hits:
        table = table.T
    table.index = table.index.map(str)
    # Keep only source rows that correspond to metadata samples.
    common = [s for s in map(str, sample_ids) if s in table.index]
    table = table.loc[common].copy()
    # Coerce feature values to numeric; non-numeric annotation columns drop out.
    table = table.apply(pd.to_numeric, errors="coerce")
    table = table.loc[:, table.notna().any(axis=0)]
    return table


def load_cohort(data_dir: str | Path, cohort_id: str):
    cdir = Path(data_dir) / "cohorts" / cohort_id
    meta = pd.read_csv(cdir / "metadata.tsv", sep="\t", low_memory=False)
    if "Sample" not in meta.columns:
        raise ValueError(f"{cohort_id}: metadata lacks a Sample column")
    meta["Sample"] = meta["Sample"].astype(str)
    taxa = orient_samples_by_features(_read_table(cdir / "genera.tsv"), meta["Sample"])
    mtb = orient_samples_by_features(_read_table(cdir / "mtb.tsv"), meta["Sample"])
    common = [s for s in meta["Sample"] if s in taxa.index and s in mtb.index]
    meta = meta.set_index("Sample").loc[common].reset_index()
    taxa = taxa.loc[common]
    mtb = mtb.loc[common]
    return taxa, mtb, meta
