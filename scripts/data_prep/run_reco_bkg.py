#!/usr/bin/env python3
"""

python3 scripts/data_prep/run_reco_bkg.py 199 299 -o data/reco_bkg_apr --list data/metadata/bkg_apr.txt --delete-intermediates

"""

import argparse
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LIST_FILE = ROOT / "data/metadata/bkg.txt"
MISSING_LOG = ROOT / "data/metadata/missing_files.txt"
DOWNLOAD_SH = ROOT / "scripts/data_prep/download.sh"
FILTER_BIN = ROOT / "build/filter"
RAW_ROOT = ROOT / "data/raw"
FILTERED_ROOT = ROOT / "data/filtered"
RECO_COLLECTIONS = (
    "EventHeader,MCParticles,LFHCALHits,LFHCALHitsContributions,"
    "LFHCALRawHits,LFHCALRawHitAssociations,LFHCALRecHits,"
    "LFHCALIslandProtoClusters,LFHCALClusters,LFHCALClusterAssociations"
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("start", type=int, help="0-based start index in the selected dataset list")
    p.add_argument("end", type=int, nargs="?", help="0-based end index in the selected dataset list")
    p.add_argument("-o", default="data/reco")
    p.add_argument("--list", default=str(DEFAULT_LIST_FILE.relative_to(ROOT)))
    p.add_argument("--delete-intermediates", action="store_true")
    return p.parse_args()


def raw_path(dataset: str) -> Path:
    scope, rel = dataset.split(":", 1)
    return RAW_ROOT / scope / rel.lstrip("/")


def tag_from_dataset(dataset: str) -> str:
    stem = Path(dataset.split(":", 1)[1]).name
    return stem.removesuffix(".edm4hep.root").split(".")[-1]


def filtered_path(dataset: str) -> Path:
    return FILTERED_ROOT / f"lfhcal_{tag_from_dataset(dataset)}.root"


def reco_path(dataset: str, output_dir: Path) -> Path:
    return output_dir / f"reco_lfhcal_{tag_from_dataset(dataset)}.root"


def run(cmd, *, check=True):
    print("+", " ".join(str(x) for x in cmd), flush=True)
    return subprocess.run(cmd, check=check, cwd=ROOT)


def append_missing(index: int, dataset: str):
    MISSING_LOG.parent.mkdir(parents=True, exist_ok=True)
    with MISSING_LOG.open("a") as out:
        out.write(f"{index}\t{dataset}\n")


def main():
    args = parse_args()
    end = args.end if args.end is not None else args.start
    output_dir = ROOT / args.o
    list_file = Path(args.list)
    if not list_file.is_absolute():
        list_file = ROOT / list_file

    datasets = list_file.read_text().splitlines()
    selected = datasets[args.start:end + 1]

    FILTERED_ROOT.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    for index, dataset in enumerate(selected, start=args.start):
        raw = raw_path(dataset)
        filtered = filtered_path(dataset)
        reco = reco_path(dataset, output_dir)

        download = run([str(DOWNLOAD_SH), str(index), str(index), str(list_file)], check=False)
        print()
        if download.returncode != 0 or not raw.exists():
            print(f"WARNING missing file at index {index}: {dataset}", flush=True)
            append_missing(index, dataset)
            continue

        run([str(FILTER_BIN), "-i", str(raw), "-o", str(filtered)])
        print()

        run([
            "eicrecon",
            f"-Ppodio:output_file={reco}",
            f"-Ppodio:output_collections={RECO_COLLECTIONS}",
            str(filtered),
        ])
        print()

        if args.delete_intermediates:
            print(f"+ rm {raw}", flush=True)
            raw.unlink()
            print(f"+ rm {filtered}", flush=True)
            filtered.unlink()
            print()


if __name__ == "__main__":
    main()
