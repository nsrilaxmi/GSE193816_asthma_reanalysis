from __future__ import annotations

import urllib.request
from pathlib import Path

from common import DATA_RAW, decompress_gzip, ensure_dirs, load_config


def download(url: str, output: Path) -> None:
    if output.exists() and output.stat().st_size > 0:
        print(f"Already present: {output.name}")
        return
    print(f"Downloading {url}")
    urllib.request.urlretrieve(url, output)


def main() -> None:
    ensure_dirs()
    config = load_config()
    base_url = config["geo_base_url"]
    for filename in config["files"].values():
        gz_name = f"{filename}.gz"
        gz_path = DATA_RAW / gz_name
        h5ad_path = DATA_RAW / filename
        if h5ad_path.exists() and h5ad_path.stat().st_size > 0:
            print(f"Already decompressed: {h5ad_path.name}")
            continue
        download(f"{base_url}/{gz_name}", gz_path)
        decompressed = decompress_gzip(gz_path, keep=True)
        print(f"Ready: {decompressed}")


if __name__ == "__main__":
    main()

