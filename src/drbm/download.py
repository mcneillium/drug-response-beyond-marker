"""Download the public data this study uses, and record checksums.

Two sources, both from the Wellcome Sanger Institute, both open with no login:

* GDSC (Genomics of Drug Sensitivity in Cancer), release 8.5: how strongly each
  drug slows the growth of each cancer cell line.
* Cell Model Passports: molecular measurements of the same cell lines (gene
  expression, protein abundance, mutations, copy number) and what each cell
  line is.

They share one identifier for a cell line (`SIDM...`), so no name matching is
needed. The original plan named the Broad Institute's DepMap portal; its
download page sits behind an interactive browser check and cannot be scripted,
so the Sanger sources are used.

Usage:
    python -m drbm.download --dest data/raw

The gene expression archive is 1.2 GB.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

FILES = {
    "GDSC2_fitted_dose_response_27Oct23.xlsx":
        "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx",
    # The older screen, used only in Experiment 7 (repeat measurements).
    "GDSC1_fitted_dose_response_27Oct23.xlsx":
        "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx",
    "model_list_latest.csv.gz": "https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz",
    "mutations_summary_20230202.zip": "https://cog.sanger.ac.uk/cmp/download/mutations_summary_20230202.zip",
    "cnv_summary_20230303.zip": "https://cog.sanger.ac.uk/cmp/download/cnv_summary_20230303.zip",
    "Proteomics_20221214.zip": "https://cog.sanger.ac.uk/cmp/download/Proteomics_20221214.zip",
    "rnaseq_merged_20250117.zip": "https://cog.sanger.ac.uk/cmp/download/rnaseq_merged_20250117.zip",
}


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, target: Path) -> None:
    """Download one file unless it is already there."""
    if target.exists() and target.stat().st_size > 0:
        return
    partial = target.with_suffix(target.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as out:
        while chunk := response.read(1 << 20):
            out.write(chunk)
    partial.replace(target)


def download(dest: Path) -> dict:
    dest.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name, url in FILES.items():
        print(f"{name} ...", flush=True)
        fetch(url, dest / name)
        manifest[name] = {"url": url, "bytes": (dest / name).stat().st_size, "sha256": sha256_of(dest / name)}
    (dest / "MANIFEST.json").write_text(json.dumps(manifest, indent=1, sort_keys=True))
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--dest", type=Path, default=Path("data/raw"))
    args = parser.parse_args(argv)
    manifest = download(args.dest)
    print(f"Done: {len(manifest)} files. Checksums in {args.dest / 'MANIFEST.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
