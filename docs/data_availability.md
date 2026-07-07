# Data Availability

This repository does not include the processed `.h5ad` files from GEO because they are large generated data files.

Download them with:

```bash
python scripts/00_download_geo.py
```

Expected files:

```text
data/raw_geo/GSE193816_all_cells_data.h5ad
data/raw_geo/GSE193816_AEC_data.h5ad
data/raw_geo/GSE193816_t_cell_data.h5ad
data/raw_geo/GSE193816_mnp_data.h5ad
data/raw_geo/GSE193816_all_data_raw_counts.h5ad
data/raw_geo/GSE193816_all_data_log_adjusted_counts.h5ad
```

Source:

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193816

Raw FASTQ files are not provided on the GEO series page for this dataset. The workflow starts from the processed supplementary `.h5ad` objects.

