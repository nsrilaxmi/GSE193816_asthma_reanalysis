#!/usr/bin/env bash
set -euo pipefail

python scripts/01_download_and_metadata.py
python scripts/02_all_cells_overview.py
python scripts/03_epithelial_signature_analysis.py
python scripts/04_tcell_analysis.py
python scripts/05_mnp_analysis.py
python scripts/06_pseudobulk_export.py

if command -v Rscript >/dev/null 2>&1; then
  Rscript scripts/06_pseudobulk_de.R
else
  echo "Rscript not found; skipping DESeq2 pseudo-bulk differential expression."
fi

if command -v quarto >/dev/null 2>&1; then
  quarto render reports/08_final_report.qmd --output-dir ../results/reports
else
  echo "Quarto not found; skipping report render."
fi

