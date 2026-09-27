# FHCal Rate study

Code for estimating hit occupancy and readout rates in the ePIC LFHCal and forward calorimeter insert from simulated EDM4hep data.

## Build

In an EIC software environment with ROOT, podio, and EDM4hep available, run:

```bash
./build.sh
```

## Run

The main programs take a ROOT file or a directory of ROOT files. For example:

```bash
./build/lfhcal_occupancy -i data/bkg_july -o lfhcal_plots/occupancy.root
./build/insert_occupancy -i data/bkg_july -o insert_plots/occupancy.root
```

The output ROOT files contain occupancy, hit-rate, and data-rate plots. The insert analysis compares the original cells with two alternative channel layouts.

Hits in the same channel are summed within each event before applying a threshold. The standard plots use 0.1 MIP for the LFHCal and 0.5 MIP for the insert; both analyses also scan other thresholds. Each accepted event represents 2 µs of simulated time, which is used to convert fired-event counts and estimated payload sizes into rates.

To download and filter datasets from a Rucio metadata list, run:

```bash
python3 scripts/data_prep/run_bkg.py 0 9 -o data/bkg_july --list data/metadata/bkg_july_minbias.txt
```

The indices select an inclusive, zero-based range from the list. This step requires a working Rucio environment. `scripts/data_prep/run_reco_bkg.py` also runs EICrecon after filtering.
