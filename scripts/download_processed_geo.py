#!/usr/bin/env python3
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

DEFAULT_URLS = {
    "gse123976": [],
    "gse197670": [],
    "gse197672": [],
    "gse116250": [],
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(DEFAULT_URLS))
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--url", action="append", default=[])
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    urls = args.url or DEFAULT_URLS[args.dataset]

    if not urls:
        print(
            "No default direct URLs are pinned. Provide one or more --url values from GEO supplementary files "
            "to download processed matrices reproducibly."
        )
        return

    for url in urls:
        target = outdir / Path(url).name
        print(f"Downloading {url} -> {target}")
        urllib.request.urlretrieve(url, target)


if __name__ == "__main__":
    main()
