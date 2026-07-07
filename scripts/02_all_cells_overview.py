from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

from common import (
    DATA_RAW,
    RESULTS_FIGURES,
    RESULTS_TABLES,
    derive_design_columns,
    ensure_dirs,
    infer_celltype_column,
    infer_sample_column,
    load_config,
)


def main() -> None:
    ensure_dirs()
    config = load_config()
    adata = sc.read_h5ad(DATA_RAW / config["files"]["all_cells"])
    sample_col = infer_sample_column(adata.obs, config)
    adata.obs = derive_design_columns(adata.obs, sample_col)
    celltype_col = infer_celltype_column(adata.obs, config)

    adata.obs["sample"] = adata.obs[sample_col].astype(str)
    adata.obs["cell_type_clean"] = adata.obs[celltype_col].astype(str)

    cell_counts = (
        adata.obs.groupby(["sample", "subject", "group", "condition"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    cell_counts.to_csv(RESULTS_TABLES / "cell_counts_by_sample.csv", index=False)

    comp = (
        adata.obs.groupby(["sample", "subject", "group", "condition", "cell_type_clean"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    totals = cell_counts[["sample", "n_cells"]].rename(columns={"n_cells": "total_cells"})
    comp = comp.merge(totals, on="sample")
    comp["fraction"] = comp["n_cells"] / comp["total_cells"]
    comp.to_csv(RESULTS_TABLES / "celltype_composition_by_sample.csv", index=False)

    if "X_umap" in adata.obsm:
        for color in ["group", "condition", "cell_type_clean"]:
            sc.pl.umap(adata, color=color, frameon=False, show=False)
            plt.savefig(RESULTS_FIGURES / f"all_cells_umap_{color}.png", dpi=300, bbox_inches="tight")
            plt.close()

    plt.figure(figsize=(10, 4))
    sns.barplot(data=cell_counts, x="sample", y="n_cells", hue="group")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig(RESULTS_FIGURES / "cell_counts_by_sample.png", dpi=300)
    plt.close()

    top_types = (
        comp.groupby("cell_type_clean", observed=True)["n_cells"]
        .sum()
        .sort_values(ascending=False)
        .head(15)
        .index
    )
    plot_comp = comp[comp["cell_type_clean"].isin(top_types)].copy()
    pivot = plot_comp.pivot_table(
        index="sample", columns="cell_type_clean", values="fraction", fill_value=0, observed=True
    )
    sample_order = cell_counts.sort_values(["group", "subject", "condition"])["sample"]
    pivot = pivot.reindex(sample_order)
    pivot.plot(kind="bar", stacked=True, figsize=(12, 5), width=0.85)
    plt.ylabel("Fraction of cells")
    plt.xlabel("")
    plt.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(RESULTS_FIGURES / "celltype_composition_stacked_bar.png", dpi=300)
    plt.close()

    print(f"Using sample column: {sample_col}")
    print(f"Using cell type column: {celltype_col}")
    print(f"Wrote overview tables and figures to {RESULTS_TABLES} and {RESULTS_FIGURES}")


if __name__ == "__main__":
    main()

