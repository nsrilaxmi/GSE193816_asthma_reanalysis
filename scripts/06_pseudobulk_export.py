from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from common import DATA_RAW, RESULTS_TABLES, derive_design_columns, ensure_dirs, infer_celltype_column, infer_sample_column, load_config


def aggregate_counts(adata, sample_col: str, celltype_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    obs = adata.obs.copy()
    obs["sample"] = obs[sample_col].astype(str)
    obs["cell_type_clean"] = obs[celltype_col].astype(str)
    obs["pseudobulk_id"] = obs["sample"] + "|" + obs["cell_type_clean"]

    groups = obs["pseudobulk_id"].astype("category")
    design = pd.get_dummies(groups, sparse=True)
    group_matrix = design.sparse.to_coo().T.tocsr()
    x = adata.X
    if not sparse.issparse(x):
        x = sparse.csr_matrix(x)
    counts = group_matrix @ x

    count_df = pd.DataFrame(
        counts.toarray(),
        index=design.columns,
        columns=adata.var_names,
    )
    meta = (
        obs[["pseudobulk_id", "sample", "subject", "group", "condition", "cell_type_clean"]]
        .drop_duplicates()
        .set_index("pseudobulk_id")
        .loc[count_df.index]
    )
    meta["n_cells"] = obs.groupby("pseudobulk_id", observed=True).size().reindex(count_df.index).to_numpy()
    return count_df, meta


def main() -> None:
    ensure_dirs()
    config = load_config()
    preferred_key = "raw_counts"
    path = DATA_RAW / config["files"][preferred_key]
    if not path.exists():
        path = DATA_RAW / config["files"]["all_cells"]
        print(f"Raw counts file not found; using {path.name}. Confirm X contains counts before DE.")

    adata = sc.read_h5ad(path)
    sample_col = infer_sample_column(adata.obs, config)
    adata.obs = derive_design_columns(adata.obs, sample_col)
    celltype_col = infer_celltype_column(adata.obs, config)

    counts, meta = aggregate_counts(adata, sample_col, celltype_col)
    counts = counts.round().astype(np.int64)
    counts.to_csv(RESULTS_TABLES / "pseudobulk_counts_sample_by_celltype.csv")
    meta.to_csv(RESULTS_TABLES / "pseudobulk_metadata_sample_by_celltype.csv")
    print(f"Exported {counts.shape[0]} pseudo-bulk profiles across {counts.shape[1]} genes.")


if __name__ == "__main__":
    main()

