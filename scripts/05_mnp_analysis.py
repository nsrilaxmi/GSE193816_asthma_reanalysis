from __future__ import annotations

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

from common import DATA_RAW, RESULTS_FIGURES, RESULTS_TABLES, derive_design_columns, ensure_dirs, infer_celltype_column, infer_sample_column, load_config


GENE_SETS = {
    "DC2_CD1C": ["CD1C", "FCER1A", "CLEC10A", "HLA-DRA", "HLA-DPB1"],
    "CCR2_MONOCYTE": ["CCR2", "S100A8", "S100A9", "FCN1", "LYZ", "CTSS"],
    "MACROPHAGE_REPAIR": ["MRC1", "MARCO", "CD163", "APOE", "C1QA", "C1QB", "MERTK"],
    "TYPE2_INFLAMMATORY_MNP": ["CCL17", "CCL22", "TNFSF13B", "IL1B", "TNF"],
}


def main() -> None:
    ensure_dirs()
    config = load_config()
    adata = sc.read_h5ad(DATA_RAW / config["files"]["mnp"])
    sample_col = infer_sample_column(adata.obs, config)
    adata.obs = derive_design_columns(adata.obs, sample_col)
    celltype_col = infer_celltype_column(adata.obs, config)
    adata.obs["sample"] = adata.obs[sample_col].astype(str)
    adata.obs["cell_type_clean"] = adata.obs[celltype_col].astype(str)

    scored = []
    for name, genes in GENE_SETS.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if len(present) >= 2:
            sc.tl.score_genes(adata, gene_list=present, score_name=name)
            scored.append(name)

    sample_scores = (
        adata.obs[["sample", "subject", "group", "condition"] + scored]
        .groupby(["sample", "subject", "group", "condition"], observed=True)
        .mean(numeric_only=True)
        .reset_index()
    )
    sample_scores.to_csv(RESULTS_TABLES / "mnp_signature_scores_by_sample.csv", index=False)

    comp = (
        adata.obs.groupby(["sample", "subject", "group", "condition", "cell_type_clean"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    totals = comp.groupby("sample", observed=True)["n_cells"].sum().rename("total_cells")
    comp = comp.merge(totals, on="sample")
    comp["fraction"] = comp["n_cells"] / comp["total_cells"]
    comp.to_csv(RESULTS_TABLES / "mnp_composition_by_sample.csv", index=False)

    if "X_umap" in adata.obsm:
        for color in ["group", "condition", "cell_type_clean"]:
            sc.pl.umap(adata, color=color, frameon=False, show=False)
            plt.savefig(RESULTS_FIGURES / f"mnp_umap_{color}.png", dpi=300, bbox_inches="tight")
            plt.close()
        genes = [gene for gene in ["CD1C", "CCR2", "FCN1", "MRC1", "CD163", "APOE"] if gene in adata.var_names]
        if genes:
            sc.pl.umap(adata, color=genes, frameon=False, show=False)
            plt.savefig(RESULTS_FIGURES / "mnp_umap_marker_genes.png", dpi=300, bbox_inches="tight")
            plt.close()

    if scored:
        long_scores = sample_scores.melt(
            id_vars=["sample", "subject", "group", "condition"],
            value_vars=scored,
            var_name="signature",
            value_name="score",
        )
        g = sns.catplot(data=long_scores, x="condition", y="score", hue="group", col="signature", kind="point", sharey=False)
        g.savefig(RESULTS_FIGURES / "mnp_signature_scores_by_condition.png", dpi=300)

    print(f"Scored MNP signatures: {', '.join(scored)}")


if __name__ == "__main__":
    main()
