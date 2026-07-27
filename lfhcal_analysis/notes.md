# LFHCAL analysis guide

`lfhcal_occupancy` reads `LFHCALHits` once and writes the channel, chip, MIP-scan,
and radius products into one ROOT file. Its structure follows `insert_analysis/occupancy`:

- `occupancy/occupancy.cc` handles files, events, decoding, and event-local channel sums.
- `occupancy/src/family.cc` applies the nominal channel threshold and routes each event.
- `occupancy/src/xy.cc` accumulates and writes origin-separated XY and YZ occupancy maps.
- `occupancy/src/channel.cc` writes channel hit-rate and data-rate products.
- `occupancy/src/chip.cc` writes chip rates, data rates, and fired-channel contours.
- `occupancy/src/mip.cc` writes channel and chip threshold scans.
- `occupancy/src/multiplicity.cc` writes hit multiplicity per passing channel.
- `occupancy/src/radius.cc` writes chip data-rate percentiles versus radius.
- `occupancy/src/central.cc` writes central-channel and origin-separated threshold products.
- `occupancy/src/shared.cc` contains shared labels and ROOT drawing helpers.

The nominal summed-channel threshold is 0.1 MIP. Threshold scans retain the full
0.0--1.5 MIP range, and representative products use 0.1, 0.3, and 0.5 MIP.

`rate_summary.cc` and `energy_spectrum.cc` remain standalone because they write a
CSV source summary and an energy spectrum instead of occupancy products.

`print_occupancy_tables.C` prints the representative threshold summary and writes
`lfhcal_plots/tables/lfhcal.csv`.

```text
./build/lfhcal_occupancy -i data/bkg_july -o lfhcal_plots/occupancy.root
```
