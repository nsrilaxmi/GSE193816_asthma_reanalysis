from __future__ import annotations

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
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


GENE_SETS = {
    "MUCUS_METAPLASIA": ["MUC5AC", "MUC5B", "SPDEF", "AGR2", "CLCA1", "POSTN"],
    "TYPE2_IL13_RESPONSE": ["POSTN", "CLCA1", "SERPINB2", "CCL26", "DPP4"],
    "MATRIX_REMODELING": ["MMP1", "MMP9", "MMP10", "MMP12", "TIMP1", "COL1A1"],
    "GLYCOLYSIS": ["SLC2A1", "HK2", "PFKP", "ALDOA", "ENO1", "LDHA"],
    "ANTIOXIDANT_RESPONSE": ["NQO1", "HMOX1", "GCLC", "GCLM", "TXNRD1", "SOD2"],
    "INJURY_REPAIR": ["AREG", "EREG", "KRT6A", "KRT16", "KRT17", "ITGA6"],
}


def main() -> None:
    ensure_dirs()
    config = load_config()
    adata = sc.read_h5ad(DATA_RAW / config["files"]["aec"])
    sample_col = infer_sample_column(adata.obs, config)
    adata.obs = derive_design_columns(adata.obs, sample_col)
    celltype_col = infer_celltype_column(adata.obs, config)
    adata.obs["sample"] = adata.obs[sample_col].astype(str)
    adata.obs["cell_type_clean"] = adata.obs[celltype_col].astype(str)

    scored = []
    for name, genes in GENE_SETS.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if len(present) < 2:
            print(f"Skipping {name}; only found {present}")
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=name)
        scored.append(name)

    scores = adata.obs[["sample", "subject", "group", "condition"] + scored].copy()
    sample_scores = scores.groupby(["sample", "subject", "group", "condition"], observed=True).mean(numeric_only=True).reset_index()
    sample_scores.to_csv(RESULTS_TABLES / "aec_signature_scores_by_sample.csv", index=False)

    if "X_umap" in adata.obsm:
        for color in ["group", "condition", "cell_type_clean"]:
            sc.pl.umap(adata, color=color, frameon=False, show=False)
            plt.savefig(RESULTS_FIGURES / f"aec_umap_{color}.png", dpi=300, bbox_inches="tight")
            plt.close()
        feature_genes = [gene for gene in ["MUC5AC", "POSTN", "CLCA1", "MMP10", "HMOX1", "AREG"] if gene in adata.var_names]
        if feature_genes:
            sc.pl.umap(adata, color=feature_genes, frameon=False, show=False)
            plt.savefig(RESULTS_FIGURES / "aec_umap_marker_genes.png", dpi=300, bbox_inches="tight")
            plt.close()

    long_scores = sample_scores.melt(
        id_vars=["sample", "subject", "group", "condition"],
        value_vars=scored,
        var_name="signature",
        value_name="score",
    )
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=long_scores, x="signature", y="score", hue="group", showfliers=False)
    sns.stripplot(data=long_scores, x="signature", y="score", hue="group", dodge=True, color="black", alpha=0.6, size=3)
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(RESULTS_FIGURES / "aec_signature_scores_by_group.png", dpi=300)
    plt.close()

    g = sns.catplot(
        data=long_scores,
        x="condition",
        y="score",
        hue="group",
        col="signature",
        kind="point",
        col_wrap=3,
        sharey=False,
        height=3,
        aspect=1.25,
    )
    g.savefig(RESULTS_FIGURES / "aec_signature_scores_by_condition.png", dpi=300)
    plt.close("all")

    print(f"Scored epithelial signatures: {', '.join(scored)}")


if __name__ == "__main__":
    main()
