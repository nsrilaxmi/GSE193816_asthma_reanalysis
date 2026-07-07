from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = PROJECT_ROOT / "data" / "raw_geo"
DATA_METADATA = PROJECT_ROOT / "data" / "metadata"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
RESULTS_FIGURES = PROJECT_ROOT / "results" / "figures"
RESULTS_TABLES = PROJECT_ROOT / "results" / "tables"


def ensure_dirs() -> None:
    for path in [DATA_RAW, DATA_METADATA, DATA_PROCESSED, RESULTS_FIGURES, RESULTS_TABLES]:
        path.mkdir(parents=True, exist_ok=True)


def load_config() -> dict:
    with (PROJECT_ROOT / "config.yaml").open() as handle:
        return yaml.safe_load(handle)


def decompress_gzip(path: Path, keep: bool = True) -> Path:
    if path.suffix != ".gz":
        return path
    output = path.with_suffix("")
    if output.exists() and output.stat().st_size > 0:
        return output
    with gzip.open(path, "rb") as src, output.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    if not keep:
        path.unlink()
    return output


def first_existing_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lower_lookup = {col.lower(): col for col in df.columns}
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
        if candidate.lower() in lower_lookup:
            return lower_lookup[candidate.lower()]
    return None


def infer_sample_column(obs: pd.DataFrame, config: dict) -> str:
    col = first_existing_column(obs, config["metadata"]["preferred_sample_columns"])
    if col:
        return col
    for candidate in obs.columns:
        values = obs[candidate].astype(str)
        if values.str.contains(r"^\d+_(AA|ANA)_(Ag|Dil|Pre)$", regex=True).mean() > 0.5:
            return candidate
    obs["_sample_from_index"] = obs.index.astype(str).str.extract(
        r"(\d+_(?:AA|ANA)_(?:Ag|Dil|Pre))", expand=False
    )
    if obs["_sample_from_index"].notna().any():
        return "_sample_from_index"
    raise ValueError("Could not infer sample column. Inspect adata.obs.columns and update config.yaml.")


def derive_design_columns(obs: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    out = obs.copy()
    sample = out[sample_col].astype(str)

    if "subject" not in out.columns:
        out["subject"] = sample.str.extract(r"(^\d+)")
    if "group" not in out.columns:
        out["group"] = sample.str.extract(r"_(AA|ANA)_")
    if "condition" not in out.columns:
        out["condition"] = sample.str.extract(r"_(Ag|Dil|Pre)$")

    out["group"] = out["group"].replace({"AA": "AA", "ANA": "ANA"})
    out["condition"] = pd.Categorical(out["condition"], categories=["Pre", "Dil", "Ag"], ordered=True)
    return out


def infer_celltype_column(obs: pd.DataFrame, config: dict) -> str:
    col = first_existing_column(obs, config["metadata"]["preferred_celltype_columns"])
    if col:
        return col
    non_design_cols = {"subject", "group", "condition"}
    categorical_cols = [
        col for col in obs.columns
        if col not in non_design_cols and (str(obs[col].dtype) == "category" or obs[col].nunique() < 80)
    ]
    if categorical_cols:
        return categorical_cols[0]
    raise ValueError("Could not infer cell type column. Inspect adata.obs.columns and update config.yaml.")


def write_sample_metadata(obs: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    cols = [sample_col, "subject", "group", "condition"]
    sample_meta = (
        obs[cols]
        .drop_duplicates()
        .rename(columns={sample_col: "sample"})
        .sort_values(["group", "subject", "condition", "sample"])
    )
    sample_meta.to_csv(DATA_METADATA / "sample_metadata_clean.csv", index=False)
    return sample_meta


def save_obs_var(adata, prefix: str) -> None:
    adata.obs.to_csv(DATA_METADATA / f"{prefix}_cell_metadata.csv")
    adata.var.to_csv(DATA_METADATA / f"{prefix}_gene_metadata.csv")

