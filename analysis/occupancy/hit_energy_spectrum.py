#!/usr/bin/env python3

"""
python3 analysis/occupancy/hit_energy_spectrum.py -i data/reco -o plots/occupancy/lfhcal_hit_energy.root
"""

import argparse
from pathlib import Path

import awkward as ak
import boost_histogram as bh
import uproot


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", required=True)
    p.add_argument("-o", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    input_dir = Path(args.i)
    files = sorted(input_dir.glob("*.root"))
    if not files:
        raise SystemExit(f"no ROOT files found in {input_dir}")

    hist = bh.Histogram(bh.axis.Regular(100, 1e-6, 1.0, transform=bh.axis.transform.log))

    for path in files:
        tree = uproot.open(f"{path}:events")
        energies = tree["LFHCALHits.energy"].array(library="ak")
        hist.fill(ak.to_numpy(ak.flatten(energies, axis=None)))

    output_path = Path(args.o)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with uproot.recreate(output_path) as fout:
        fout["lfhcal_hit_energy"] = hist


if __name__ == "__main__":
    main()
