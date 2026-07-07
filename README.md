# Reproducible GEO Reanalysis of Allergen-Induced Airway Remodeling in Asthma

This repository contains a reproducible reanalysis workflow for GEO dataset `GSE193816`, a human single-cell RNA-seq study of endobronchial brush samples from allergic asthmatics and allergic non-asthmatic controls before and after airway challenge.

The goal is to provide a clean, independent, GitHub-ready analysis that can be run from processed GEO `.h5ad` files. It focuses on sample metadata curation, cell composition, airway epithelial remodeling signatures, T-cell and mononuclear phagocyte programs, and pseudo-bulk differential expression inputs.

## Biological Question

How does allergen exposure change airway epithelial, T-cell, and mononuclear phagocyte programs differently in allergic asthmatics compared with allergic non-asthmatic controls?

Main hypothesis tested:

> In response to allergen, asthmatic airways show stronger epithelial remodeling and type-2 inflammatory programs, while allergic non-asthmatic controls show more repair/pro-resolution responses.

Published findings are treated as hypotheses to reproduce and evaluate, not as assumed results.

## Repository Status

This repository is intentionally safe for public GitHub upload:

- Large GEO `.h5ad` data files are ignored by Git.
- Generated figures, tables, and rendered reports are ignored by default.
- Placeholder `.gitkeep` files preserve the expected output folders.
- Download, analysis, and report scripts are tracked.
- Lightweight GitHub Actions checks validate Python syntax and R script parsing.

## Project Layout

```text
GSE193816_asthma_reanalysis/
├── .github/workflows/
├── data/
│   ├── raw_geo/
│   ├── processed/
│   └── metadata/
├── docs/
├── scripts/
├── notebooks/
├── results/
│   ├── figures/
│   ├── tables/
│   └── reports/
├── reports/
├── config.yaml
├── environment.yml
├── requirements.txt
└── README.md
```

## Data

GEO accession: `GSE193816`

The workflow uses processed supplementary `.h5ad` files from GEO. Raw FASTQ files are not included in GEO for this series, and the processed data are not committed to this repository because they are large.

The main files expected in `data/raw_geo/` are:

```text
GSE193816_all_cells_data.h5ad
GSE193816_AEC_data.h5ad
GSE193816_t_cell_data.h5ad
GSE193816_mnp_data.h5ad
GSE193816_all_data_raw_counts.h5ad
GSE193816_all_data_log_adjusted_counts.h5ad
```

## Installation

```bash
conda env create -f environment.yml
conda activate gse193816
```

Or with pip:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Reproduce the Analysis

From the repository root:

```bash
python scripts/00_download_geo.py
bash scripts/run_all.sh
```

Or run each step manually:

```bash
python scripts/01_download_and_metadata.py
python scripts/02_all_cells_overview.py
python scripts/03_epithelial_signature_analysis.py
python scripts/04_tcell_analysis.py
python scripts/05_mnp_analysis.py
python scripts/06_pseudobulk_export.py
```

If R/DESeq2 is available:

```bash
Rscript scripts/06_pseudobulk_de.R
```

Render the report:

```bash
quarto render reports/08_final_report.qmd --output-dir ../results/reports
```

## Primary Outputs

The analysis creates the following local outputs:

```text
data/metadata/
├── gse193816_all_cells_cell_metadata.csv
├── gse193816_all_cells_gene_metadata.csv
└── sample_metadata_clean.csv

results/tables/
├── cell_counts_by_sample.csv
├── celltype_composition_by_sample.csv
├── aec_signature_scores_by_sample.csv
├── tcell_signature_scores_by_sample.csv
├── mnp_signature_scores_by_sample.csv
├── pseudobulk_counts_sample_by_celltype.csv
└── pseudobulk_metadata_sample_by_celltype.csv

results/figures/
├── all_cells_umap_group.png
├── all_cells_umap_condition.png
├── celltype_composition_stacked_bar.png
├── aec_signature_scores_by_condition.png
└── additional lineage-specific UMAP/signature plots
```

Generated outputs are ignored by Git so the repository stays small. For a portfolio release, selected final figures can be copied into `docs/figures/` and explicitly committed.

## Analysis Modules

- `01_download_and_metadata.py`: loads all-cell object, infers sample/group/condition fields, exports metadata.
- `02_all_cells_overview.py`: UMAP overview, cell counts, and composition tables.
- `03_epithelial_signature_analysis.py`: airway epithelial remodeling signatures and marker feature plots.
- `04_tcell_analysis.py`: TH2/IL9 and activation signature summaries.
- `05_mnp_analysis.py`: DC2, CCR2 monocyte, macrophage repair, and inflammatory MNP summaries.
- `06_pseudobulk_export.py`: exports sample-by-cell-type pseudo-bulk count matrices.
- `06_pseudobulk_de.R`: DESeq2 template for patient-aware contrasts where replicate structure allows.

## Limitations

The cohort is small and GEO lists 21 samples, so not every subject has every condition. Composition and signature results should be interpreted as exploratory. Pseudo-bulk contrasts should be checked for sufficient biological replication before strong biological claims are made.

## Author

**Srilaxmi Nerella**  
UCSF profile: [profiles.ucsf.edu/srilaxmi.nerella](https://profiles.ucsf.edu/srilaxmi.nerella)  
Google Scholar: [scholar.google.com/citations?user=wjN338cAAAAJ](https://scholar.google.com/citations?user=wjN338cAAAAJ&hl=en)

## References

- GEO accession: [GSE193816](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193816)
- Original project resources: [Villani Lab asthma browser](https://villani.mgh.harvard.edu/allergy-asthma/)
- Original code resource: [villani-lab/airway_allergic_asthma](https://github.com/villani-lab/airway_allergic_asthma)
