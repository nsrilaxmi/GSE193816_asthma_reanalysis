from __future__ import annotations

import scanpy as sc

from common import DATA_RAW, derive_design_columns, ensure_dirs, infer_sample_column, load_config, save_obs_var, write_sample_metadata


def main() -> None:
    ensure_dirs()
    config = load_config()
    path = DATA_RAW / config["files"]["all_cells"]
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}. Run scripts/00_download_geo.py first.")

    adata = sc.read_h5ad(path)
    sample_col = infer_sample_column(adata.obs, config)
    adata.obs = derive_design_columns(adata.obs, sample_col)

    save_obs_var(adata, "gse193816_all_cells")
    sample_meta = write_sample_metadata(adata.obs, sample_col)

    print(adata)
    print("obs columns:", adata.obs.columns.tolist())
    print("obsm keys:", list(adata.obsm.keys()))
    print(sample_meta)


if __name__ == "__main__":
    main()

