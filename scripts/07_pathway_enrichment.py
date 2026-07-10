from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from common import PROJECT_ROOT, RESULTS_FIGURES, RESULTS_TABLES, ensure_dirs


DOCS_FIGURES = PROJECT_ROOT / "docs" / "figures"
DOCS_TABLES = PROJECT_ROOT / "docs" / "tables"


GENE_SETS = {
    "Type 2 epithelial response": {
        "POSTN", "CLCA1", "SERPINB2", "SERPINB4", "CCL26", "DPP4", "SLC26A4", "ALOX15"
    },
    "Mucus/metaplasia": {
        "MUC5AC", "MUC5B", "MUC2", "SPDEF", "AGR2", "FOXA3", "CLCA1", "TFF3"
    },
    "Matrix remodeling": {
        "TNC", "MMP1", "MMP7", "MMP9", "MMP10", "MMP12", "TIMP1", "COL1A1", "COL3A1", "POSTN"
    },
    "Epithelial injury/repair": {
        "AREG", "EREG", "KRT6A", "KRT16", "KRT17", "ITGA6", "CD55", "RARRES1", "TNC"
    },
    "Ciliated epithelial program": {
        "FOXJ1", "PIFO", "TPPP3", "SNTN", "DNAH5", "DNAH9", "CAPS", "RSPH1"
    },
    "Oxidative stress response": {
        "HMOX1", "NQO1", "GCLC", "GCLM", "TXNRD1", "SOD2", "FTH1", "FTL"
    },
    "T-cell activation": {
        "CD69", "IL2RA", "ICOS", "TNFRSF4", "CD40LG", "NFKBIA", "JUN", "FOS"
    },
    "Type 2 T-cell program": {
        "GATA3", "IL4", "IL5", "IL9", "IL13", "CCR4", "PTGDR2", "ICOS"
    },
    "Antigen presentation": {
        "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA", "CD86", "CD40"
    },
    "Monocyte/DC inflammatory program": {
        "FCN1", "CCR2", "S100A8", "S100A9", "LYZ", "IL1B", "TNF", "CCL2", "CCL3", "CCL4"
    },
    "Macrophage lipid/complement program": {
        "APOC1", "APOE", "SCD", "C1QA", "C1QB", "C1QC", "MARCO", "MRC1", "CD163"
    },
}


CONTRAST_LABELS = {
    "pseudobulk_AEC_AA_Ag_vs_Pre_top50.csv": "AEC AA Ag vs Pre",
    "pseudobulk_AEC_ANA_Ag_vs_Pre_top50.csv": "AEC ANA Ag vs Pre",
    "pseudobulk_Tcell_AA_Ag_vs_Pre_top50.csv": "T cell AA Ag vs Pre",
    "pseudobulk_MNP_AA_Ag_vs_Pre_top50.csv": "MNP AA Ag vs Pre",
}


def read_contrast_table(filename: str) -> pd.DataFrame:
    for folder in [RESULTS_TABLES, DOCS_TABLES]:
        path = folder / filename
        if path.exists():
            df = pd.read_csv(path)
            df["source_file"] = filename
            return df
    raise FileNotFoundError(f"Could not find {filename} in results/tables or docs/tables")


def benjamini_hochberg(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["padj"] = multipletests(out["pvalue"], method="fdr_bh")[1]
    return out


def enrich_gene_list(
    genes: set[str],
    background: set[str],
    contrast: str,
    direction: str,
) -> list[dict]:
    rows = []
    background_size = len(background)
    query_size = len(genes)
    for pathway, pathway_genes in GENE_SETS.items():
        pathway_in_background = pathway_genes & background
        overlap = genes & pathway_genes
        if not overlap:
            continue
        pvalue = hypergeom.sf(
            len(overlap) - 1,
            background_size,
            len(pathway_in_background),
            query_size,
        )
        rows.append(
            {
                "contrast": contrast,
                "direction": direction,
                "pathway": pathway,
                "n_query_genes": query_size,
                "n_pathway_genes_in_background": len(pathway_in_background),
                "n_overlap": len(overlap),
                "overlap_genes": ";".join(sorted(overlap)),
                "overlap_fraction": len(overlap) / max(query_size, 1),
                "pvalue": pvalue,
            }
        )
    return rows


def plot_enrichment(enrichment: pd.DataFrame) -> None:
    if enrichment.empty:
        print("No pathway overlaps to plot.")
        return
    plot_df = (
        enrichment.sort_values(["n_overlap", "overlap_fraction", "pvalue"], ascending=[False, False, True])
        .groupby(["contrast", "direction"], as_index=False)
        .head(4)
        .copy()
    )
    plot_df["contrast_direction"] = plot_df["contrast"] + " | " + plot_df["direction"]
    plot_df["minus_log10_p"] = -plot_df["pvalue"].clip(lower=1e-300).map(lambda x: __import__("math").log10(x))

    plt.figure(figsize=(11, max(4.8, 0.42 * len(plot_df))))
    ax = sns.scatterplot(
        data=plot_df,
        x="minus_log10_p",
        y="pathway",
        size="n_overlap",
        hue="contrast_direction",
        sizes=(60, 260),
        alpha=0.85,
    )
    ax.set_xlabel("-log10 exploratory overlap p-value")
    ax.set_ylabel("")
    ax.legend(title="Contrast / direction", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    plt.tight_layout()
    for folder in [RESULTS_FIGURES, DOCS_FIGURES]:
        folder.mkdir(parents=True, exist_ok=True)
        plt.savefig(folder / "pathway_enrichment_dotplot.png", dpi=300, bbox_inches="tight")
    plt.close()


def write_summary(enrichment: pd.DataFrame) -> None:
    lines = [
        "# Pathway Enrichment Summary",
        "",
        "This first-pass enrichment layer tests whether committed pseudo-bulk top-gene lists overlap curated airway, immune, epithelial, and myeloid gene modules.",
        "",
        "Because the committed inputs are compact top-50 preview tables rather than full ranked differential-expression results, p-values should be treated as exploratory overlap statistics rather than formal genome-wide enrichment claims.",
        "",
        "## Top Overlaps",
        "",
    ]
    if enrichment.empty:
        lines.append("No curated pathway overlaps were detected in the available top-gene tables.")
    else:
        top = enrichment.sort_values(["n_overlap", "overlap_fraction", "pvalue"], ascending=[False, False, True]).head(12)
        for row in top.itertuples(index=False):
            gene_word = "gene" if row.n_overlap == 1 else "genes"
            lines.append(
                f"- **{row.contrast} ({row.direction})**: {row.pathway} "
                f"overlaps {row.n_overlap} {gene_word} (`{row.overlap_genes}`)."
            )
    lines.extend(
        [
            "",
            "## Recommended Interpretation",
            "",
            "- Use these results to prioritize biological themes for follow-up.",
            "- Re-run enrichment on full ranked pseudo-bulk DE results before making strong pathway-level claims.",
            "- Interpret broad-lineage T-cell and MNP results cautiously because rare subtype programs can be diluted.",
            "",
        ]
    )
    (PROJECT_ROOT / "docs" / "pathway_enrichment_summary.md").write_text("\n".join(lines))


def main() -> None:
    ensure_dirs()
    DOCS_TABLES.mkdir(parents=True, exist_ok=True)
    DOCS_FIGURES.mkdir(parents=True, exist_ok=True)

    contrast_tables = {filename: read_contrast_table(filename) for filename in CONTRAST_LABELS}
    observed_genes = set()
    for table in contrast_tables.values():
        observed_genes.update(table["gene"].dropna().astype(str).str.upper())
    pathway_background = set().union(*GENE_SETS.values())
    background = observed_genes | pathway_background

    rows = []
    for filename, table in contrast_tables.items():
        contrast = CONTRAST_LABELS[filename]
        table = table.copy()
        table["gene"] = table["gene"].astype(str).str.upper()
        table = table.dropna(subset=["log2FoldChange"])
        up_genes = set(table.loc[table["log2FoldChange"] > 0, "gene"])
        down_genes = set(table.loc[table["log2FoldChange"] < 0, "gene"])
        rows.extend(enrich_gene_list(up_genes, background, contrast, "up in allergen"))
        rows.extend(enrich_gene_list(down_genes, background, contrast, "down in allergen"))

    enrichment = benjamini_hochberg(pd.DataFrame(rows))
    if not enrichment.empty:
        enrichment = enrichment.sort_values(["padj", "pvalue", "contrast", "direction", "pathway"])

    for folder in [RESULTS_TABLES, DOCS_TABLES]:
        folder.mkdir(parents=True, exist_ok=True)
        enrichment.to_csv(folder / "pathway_enrichment_results.csv", index=False)

    plot_enrichment(enrichment)
    write_summary(enrichment)
    print(f"Wrote {len(enrichment)} pathway overlap rows.")


if __name__ == "__main__":
    main()
