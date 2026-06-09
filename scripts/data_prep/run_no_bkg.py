#!/usr/bin/env python3
"""

python3 scripts/data_prep/run_no_bkg.py 0 99 -o data/reco_no_bkg_apr --delete-intermediates

"""

import argparse
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LIST_FILE = ROOT / "data/metadata/no_bkg_apr.txt"
MISSING_LOG = ROOT / "data/metadata/missing_files.txt"
DOWNLOAD_SH = ROOT / "scripts/data_prep/download.sh"
FILTER_BIN = ROOT / "build/filter"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("start", type=int, help="0-based start index in the selected dataset list")
    p.add_argument("end", type=int, nargs="?", help="0-based end index in the selected dataset list")
    p.add_argument("-o", default="data/reco_no_bkg")
    p.add_argument("--list", default=str(DEFAULT_LIST_FILE.relative_to(ROOT)))
    p.add_argument("--delete-intermediates", action="store_true")
    return p.parse_args()


def raw_path(dataset: str) -> Path:
    scope, rel = dataset.split(":", 1)
    return ROOT / "data/raw" / scope / rel.lstrip("/")


def output_path(dataset: str, output_dir: Path) -> Path:
    return output_dir / Path(dataset.split(":", 1)[1]).name


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

    output_dir.mkdir(parents=True, exist_ok=True)

    for index, dataset in enumerate(selected, start=args.start):
        raw = raw_path(dataset)
        output = output_path(dataset, output_dir)

        download = run([str(DOWNLOAD_SH), str(index), str(index), str(list_file)], check=False)
        print()
        if download.returncode != 0 or not raw.exists():
            print(f"WARNING missing file at index {index}: {dataset}", flush=True)
            append_missing(index, dataset)
            continue

        run([str(FILTER_BIN), "-i", str(raw), "-o", str(output)])
        print()

        if args.delete_intermediates:
            print(f"+ rm {raw}", flush=True)
            raw.unlink()
            print()


if __name__ == "__main__":
    main()
