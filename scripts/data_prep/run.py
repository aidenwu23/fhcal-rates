#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LIST_FILE = ROOT / "data/metadata/file_list.txt"
DOWNLOAD_SH = ROOT / "scripts/data_prep/download.sh"
FILTER_BIN = ROOT / "build/filter"
RAW_ROOT = ROOT / "data/raw"
FILTERED_ROOT = ROOT / "data/filtered"
RECO_COLLECTIONS = (
    "EventHeader,MCParticles,LFHCALHits,LFHCALHitsContributions,"
    "LFHCALRawHits,LFHCALRawHitAssociations,LFHCALRecHits,"
    "LFHCALIslandProtoClusters,LFHCALClusters"
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("start", type=int, help="0-based start index in file_list.txt")
    p.add_argument("end", type=int, nargs="?", help="0-based end index in file_list.txt")
    p.add_argument("-o", default="data/reco")
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


def run(cmd):
    print("+", " ".join(str(x) for x in cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=ROOT)


def main():
    args = parse_args()
    end = args.end if args.end is not None else args.start
    output_dir = ROOT / args.o

    datasets = Path(LIST_FILE).read_text().splitlines()
    selected = datasets[args.start:end + 1]

    FILTERED_ROOT.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    for index, dataset in enumerate(selected, start=args.start):
        raw = raw_path(dataset)
        filtered = filtered_path(dataset)
        reco = reco_path(dataset, output_dir)

        run([str(DOWNLOAD_SH), str(index), str(index)])
        print()

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
