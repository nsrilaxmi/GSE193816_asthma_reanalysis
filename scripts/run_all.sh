#!/usr/bin/env bash
set -euo pipefail

export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/gse193816_mpl}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp/gse193816_numba}"

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_BIN" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
  elif command -v python >/dev/null 2>&1; then
    PYTHON_BIN="python"
  else
    echo "No Python interpreter found. Install Python or activate the conda environment first." >&2
    exit 1
  fi
fi

"$PYTHON_BIN" scripts/01_download_and_metadata.py
"$PYTHON_BIN" scripts/02_all_cells_overview.py
"$PYTHON_BIN" scripts/03_epithelial_signature_analysis.py
"$PYTHON_BIN" scripts/04_tcell_analysis.py
"$PYTHON_BIN" scripts/05_mnp_analysis.py
"$PYTHON_BIN" scripts/06_pseudobulk_export.py

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
