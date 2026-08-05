# LFHCAL and Forward-Insert Rate Methodology

## End-to-end flow

```text
for each requested campaign file:
    download the EDM4hep ROOT file
    retain EventHeader, MCParticles, LFHCAL collections, and insert collections
    optionally run eicrecon on the filtered file

for each analysis executable:
    find input ROOT files and sort them by path
    for each readable event containing the requested hit collection:
        increment accepted_event_count
        decode hit cell IDs
        group hits into event-local channels
        sum channel energy before applying thresholds
        update full-sample channel, chip, spatial, and source accumulators

    simulated_time = accepted_event_count * 2 microseconds
    convert accumulated event counts and payload bits into rates
    write histograms, graphs, canvases, or CSV summaries
```

The rate denominator therefore represents consecutive 2 microsecond sample windows:

```text
total_time_s = accepted_event_count * 2e-6
hit_rate_Hz = passing_event_count / total_time_s
data_rate_Gbps = accumulated_payload_bits / total_time_s / 1e9
```

An event enters the denominator when its frame is readable and contains the analyzed collection. Empty events containing that collection still enter the denominator.

## Data preparation

### Dataset runners

`run_bkg.py` and `run_no_bkg.py` perform this sequence for an inclusive, zero-based range of lines in a metadata file:

```text
read dataset names from metadata[start : end + 1]
for each dataset:
    derive its data/raw path from "scope:path"
    call download.sh for that metadata line
    if the download failed or the expected file is absent:
        append index and dataset name to data/metadata/missing_files.txt
        continue
    run build/filter on the raw file
    write the slimmed file to the selected output directory
    delete the raw file when --delete-intermediates was supplied
```

`run_reco_bkg.py` adds reconstruction:

```text
download raw file
filter raw file to data/filtered/lfhcal_<tag>.root
run eicrecon on the filtered file
write EventHeader, MCParticles, LFHCALHits, LFHCAL hit contributions,
      raw hits, associations, reconstructed hits, proto-clusters, clusters,
      and cluster associations
optionally delete raw and filtered intermediates
```

`download.sh` converts the requested zero-based indices to one-based metadata line numbers and calls `rucio download`. With `--output`, it flattens Rucio's nested output path into the requested directory.

### Collection filter

`filter` copies every event while retaining collections whose names satisfy:

```text
name == EventHeader
or name == MCParticles
or name contains "LFHCAL"
or name contains "HcalEndcapPInsert"
```

## Shared decoding and source classification

### LFHCAL cell identity

The LFHCAL decoder reads this packed bit layout:

```text
system:8, moduleIDx:6, moduleIDy:6, moduletype:1, passive:1,
towerx:2, towery:1, rlayerz:4, layerz:4
```

An analysis channel is identified by:

```text
(moduleIDx, moduleIDy, moduletype, towerx, towery, rlayerz)
```

`rlayerz` is the seven-bin longitudinal readout layer. `layerz` is decoded by the utility but excluded from analysis-channel identity. Passive cells are rejected by `lfhcal_occupancy`.

The transverse channel center is:

```text
8M module, moduletype == 0:
    x_cm = 270 - 10*moduleIDx + 7.5 - 5*towerx

4M module, moduletype == 1:
    x_cm = 265 - 10*moduleIDx + 2.5 - 5*towerx

both module types:
    y_cm = 265 - 10*moduleIDy + 2.5 - 5*towery

x_mm = 10*x_cm
y_mm = 10*y_cm
```

An LFHCAL hardware chip is identified by:

```text
(moduleIDx, moduleIDy, moduletype, module_half)

for an 8M module:
    module_half = 0 for towerx 0 or 1
    module_half = 1 for towerx 2 or 3

for a 4M module:
    module_half = 0
```

The chip identity excludes `rlayerz`, so passing channels from multiple readout layers can contribute to one chip payload in an event.

### Insert cell identity

The insert decoder reads:

```text
system:8, side:1, layer:8, slice:7, unused bits through bit 31,
x:16 signed, y:16 signed
```

The occupancy analysis accepts insert layers 1 through 60 and sides 0 and 1. It uses the simulated hit position for transverse mapping.

### Generator-status families

Every positive-energy hit contribution is assigned from its particle's `generatorStatus`:

| Status range      | Origin index | Label                   |
| ---               | ---:         | ---                     |
| `[0, 1000)`       | 0            | DIS                     |
| `[2000, 3000)`    | 1            | synchrotron radiation   |
| `[3000, 4000)`    | 2            | electron bremsstrahlung |
| `[4000, 5000)`    | 3            | electron Touschek       |
| `[5000, 6000)`    | 4            | electron Coulomb        |
| `[6000, 7000)`    | 5            | proton beam gas         |
| every other value | 6            | other/undocumented      |

Status values at the start of a shifted block, such as 2000 or 6000, remain in that broad family, since they are shower secondaries and contribute to energy depositions.

The four-class summary used by `lfhcal_rate_summary` is:

```text
[0, 1000)       -> DIS
[2000, 6000)    -> electron-beam background
[6000, 7000)    -> proton-beam background
everything else -> other
```

## LFHCAL occupancy analysis

Main executable:

```text
./build/lfhcal_occupancy -i INPUT.root|DIRECTORY -o OUTPUT.root
```

### Event construction

`lfhcal_occupancy` reads the collection named `LFHCALHits`.

```text
for each accepted event:
    event_channels = map keyed by decoded LFHCAL channel

    for each LFHCAL simulated hit:
        decode cell ID
        reject passive cells
        reject rlayerz outside [0, 6]

        channel = event_channels[channel_id]
        on the channel's first hit:
            store channel ID, chip ID, and decoded x/y center

        channel.energy += hit.energy
        channel.hit_count += 1

        for each positive-energy contribution:
            origin = family(contribution.particle.generatorStatus)
            channel.energy_by_origin[origin] += contribution.energy

    pass the completed channel list to every LFHCAL analysis product
```

Thresholds always act on the total event-local channel energy. Individual simulated hits and individual contributions are never thresholded in this executable.

### MIP definitions and channel types

```text
readout layers 0 and 1:
    channel type = 5_tile
    1 MIP = 3.5e-3 GeV

readout layers 2 through 6:
    channel type = 10_tile
    1 MIP = 7.0e-3 GeV
```

The fixed-threshold products use `0.1 MIP`. The scans use coefficients `0.0, 0.1, ..., 1.5` and the following comparison:

```text
channel passes when summed_energy_GeV > coefficient * layer_MIP_GeV
```

### Fixed-threshold channel products

```text
passing_channels = channels with energy > 0.1 * layer_MIP

for each passing channel:
    channel.passing_events += 1
    channel.total_hits += event_channel.hit_count
    increment the event's fired-channel count for its channel type
```

Outputs include fired channels per event, simulated-hit multiplicity per passing channel, and channel-rate distributions. The channel data-rate estimate assigns each passing channel:

```text
channel_payload_per_passing_event = 32 bits/hit * 4 samples = 128 bits
channel_data_rate = passing_events * 128 / total_time_s / 1e9
```

The channel-level data rate excludes the shared 128-bit chip framing term.

### Spatial products

For each fixed-threshold passing channel, the code updates its physical readout layer and the inclusive `sum_layers` group.

```text
inclusive spatial group:
    always increment

origin spatial group:
    increment when channel.energy_by_origin[origin] > 0
```

Origin maps indicate source presence inside a channel event selected by total channel energy. A channel event containing several origins increments several origin maps.

The XY maps report:

```text
occupancy = passing_channel_events / accepted_events
rate_Hz = passing_channel_events / total_time_s
```

The inclusive layer map combines all readout layers at shared transverse coordinates. YZ products use channels with `-400 mm <= x <= 300 mm`, then aggregate by readout layer and y coordinate. Their values are normalized by event count and by the number of contributing channels in each aggregate bin.

### Fixed-threshold chip products

```text
for each event:
    count passing channels assigned to each LFHCAL chip

    for each fired chip:
        chip.passing_events += 1
        chip.total_active_channels += active_channel_count
        chip.payload_bits += (128 + 32*active_channel_count) * 4
```

Outputs report chip hit rate, chip data rate, chips per event, active channels per chip-event, and the relationship between chip rate and mean fired-channel count.

### Threshold scans

For each coefficient from 0.0 through 1.5:

```text
for each event channel passing this coefficient:
    channel_passes[channel] += 1
    active_channels[chip] += 1

for each fired chip:
    chip_passes[chip] += 1
    chip_payload_bits[chip] += (128 + 32*active_channels[chip]) * 4
```

The resulting graphs contain:

- p95 and p99 channel hit rate;
- p95 and p99 chip hit rate;
- p95 and p99 chip data rate;
- mean chip data rate across chips present in the accumulated map;
- total detector data rate from the sum of all chip payloads;
- mean active channels per unthresholded chip instance and its equivalent fired-channel rate.

At each coefficient, percentile populations contain channels or chips that passed that coefficient at least once. Zero-rate geometry entries are absent. Percentile index selection is `floor(percentile * (N - 1))` after sorting.

For the mean-active-channel curve, the numerator counts threshold-passing channels while the denominator counts chip instances containing any event channel before that threshold is applied.

### Radial chip study

At `0.1`, `0.3`, and `0.5 MIP`:

```text
build framed payload per chip as above
compute chip radius from its decoded center and detector origin
place chip in one of 30 uniform bins over 0 to 3000 mm
within each populated radius bin, compute mean, p95, and p99 data rate
```

### Central-region study

The central selection is `sqrt(x_mm^2 + y_mm^2) <= 1000 mm`.

```text
for each central event channel and every MIP coefficient:
    if total channel energy passes:
        increment all-source channel passes
        increment DIS passes when DIS contributed positive energy
        increment pBeamGas passes when proton beam gas contributed positive energy
```

It writes p95/p99 threshold scans for both channel types and source groups `all`, `DIS`, and `pBeamGas`, plus full rate distributions at `0.1`, `0.3`, and `0.5 MIP`.

## Forward-insert occupancy analysis

Main executable:

```text
./build/insert_occupancy -i INPUT.root|DIRECTORY -o OUTPUT.root
```

### Event construction shared by all layouts

The executable chooses the first simulated-hit collection whose name contains `HcalEndcapPInsert` and excludes names containing `Contributions`.

```text
for each accepted event:
    event_hits = []
    for each insert hit:
        decode cell ID
        require physical layer 1 through 60
        require side 0 or 1
        append {cell_id, layer, side, hit.x_mm, hit.y_mm, hit.energy_GeV}

    send the same event_hits to six accumulators:
        three channel layouts
        times two inner-ring treatments
```

### Layout A: original insert

```text
channel identity = simulated insert cell ID
channel layer = physical_layer - 1
channel center = simulated hit position
1 MIP per channel = 4e-4 GeV
```

All 60 physical layers remain separate. Energy from hits sharing the same simulated cell is summed inside each event.

### Layout B: LFHCAL-sized tiles

```text
ix = floor(hit.x_mm / 50 mm)
iy = floor(hit.y_mm / 50 mm)
channel identity = (physical_layer - 1, ix, iy)
channel center = ((ix + 0.5)*50 mm, (iy + 0.5)*50 mm)
1 MIP per channel = 4e-4 GeV
```

All 60 physical layers remain separate. Hits inside the same 50 mm by 50 mm virtual tile and physical layer are summed inside each event.

### Layout C: full-LFHCAL segmentation

```text
physical layers  1-5  -> segment 0
physical layers  6-10 -> segment 1
physical layers 11-20 -> segment 2
physical layers 21-30 -> segment 3
physical layers 31-40 -> segment 4
physical layers 41-50 -> segment 5
physical layers 51-60 -> segment 6

ix = floor(hit.x_mm / 50 mm)
iy = floor(hit.y_mm / 50 mm)
channel identity = (segment, ix, iy)
channel_MIP = number_of_physical_layers_in_segment * 4e-4 GeV
```

Hits at one transverse location are summed across every physical layer in their segment before thresholding.

### Inner-ring treatments

Every layout is analyzed both with the complete channel set and with an inner-ring removal. The removal uses the channel center:

```text
dx = channel_center_x_mm - (-72 mm)
keep channel when dx^2 + channel_center_y_mm^2 >= (163 mm)^2
```

Original-insert channels use their simulated centers. Both virtual-tile layouts use their 50 mm tile centers. This mask is applied before occupancy, threshold scans, chip assignment, radial summaries, and side totals.

### Insert chip mappings

The original-insert and LFHCAL-tile layouts use angular chips. Chip identity is `(layer, side, angular_region)`.

```text
angle about x = -12.725 mm
convert angle to a side-local interval from 0 to pi
side 0 -> four angular regions
side 1 -> three angular regions
preserve layer in chip identity
```

The implemented region boundaries in radians are:

```text
side 0: 0, 0.9127426519, pi/2, 2.2288500017, pi
side 1: 0, 0.8133195162, 2.3282731374, pi
```

The full-LFHCAL layout groups virtual tile indices in 2 by 2 blocks:

```text
chip_ix = ix / 2 for ix >= 0, else (ix - 1) / 2
chip_iy = iy / 2 for iy >= 0, else (iy - 1) / 2
chip identity = (1, chip_ix, chip_iy)
```

The first chip-identity field is fixed to `1` by the current mapper, so the seven longitudinal segments share a transverse virtual chip identity. Each virtual chip covers 100 mm by 100 mm transversely.

### Insert thresholds and fixed products

Each layout keeps event-local hit multiplicity and energy separately:

```text
channel_hit_count[channel] += 1 per simulated hit
channel_energy[channel] += hit.energy

fixed-product channel passes when:
    channel_energy > 0.5 * channel_MIP
```

Passing-event counts determine channel hit rates. `total_hits` retains the number of simulated hits summed into those passing channel events.

Every layout writes per-layer or per-segment occupancy, an inclusive group, channel hit-rate and channel data-rate distributions, chip products, MIP scans, side data-rate scans, and chip fired-channel contours. The full-LFHCAL layout also writes radial chip summaries.

### Insert MIP scans and payload model

All layouts scan coefficients `0.0, 0.1, ..., 1.5` with strict `energy > coefficient * channel_MIP` selection.

```text
channel payload for one passing channel-event:
    32 bits * 4 samples = 128 bits

chip payload for one fired chip-event:
    (128 framing bits + 32 bits * active_channel_count) * 4 samples
```

For each coefficient, the code writes p95 and p99 channel hit rates, channel data rates, chip hit rates, and chip data rates. Percentile populations contain IDs that pass at least once, and indexing is `floor(percentile * (N - 1))` after sorting.

Side totals rebuild chip activity separately for side 0 and side 1, sum all framed chip payload in every event, and report total Gbps versus MIP coefficient.

The full-LFHCAL radial study uses coefficients `0.1`, `0.5`, and `1.5 MIP`, eight uniform radius bins from 0 to 750 mm, and p95/p99 chip data rates in each populated bin. Radius is measured from the detector origin using the 100 mm virtual-chip center.

### Insert ROOT directory structure

```text
original_insert/
    with_inner_ring/
    no_inner_ring/
lfhcal_tiles/
    with_inner_ring/
    no_inner_ring/
full_lfhcal/
    with_inner_ring/
    no_inner_ring/
```

## Standalone analyses

### LFHCAL rate summary

`lfhcal_rate_summary` groups active LFHCAL hits by event channel, sums hit energy, and applies `0.1 MIP` using `3.5e-3 GeV` for readout layers 0-1 and `7e-3 GeV` for layers 2-6.

```text
for each passing channel-event:
    increment all_sources
    inspect every positive contribution in its component hits
    increment DIS when any DIS contribution is present
    increment proton_beam_background when any pBeamGas contribution is present
    increment other_backgrounds when any remaining contribution is present
```

The source columns are presence categories and can overlap. The CSV contains p95 and p99 channel rates for all sources and the three source groups.

### LFHCAL hit origins

`hit_origins` works at simulated-hit contribution level:

```text
for each hit:
    sum positive contribution energy separately for all seven status families
    for each family with positive summed energy:
        fill that family's energy-deposit spectrum for the decoded readout layer
        fill inclusive all-layer spectrum
        fill hit-position eta when family energy passes 0 or 5e-4 GeV
        increment that family's hit count
```

It also records the raw generator-status distribution and writes an origin CSV with hit counts, fractions, and rates normalized by the 2 microsecond event window. A multi-origin hit contributes once to each represented family.

`MC_eta_and_energy` classifies every `MCParticles` entry by generator status and fills particle energy and momentum-pseudorapidity spectra.

`multiplicity` records, per LFHCAL simulated hit, the number of unique contributing MC particle object IDs and the number of unique generator-status families.

`pbeamgas_unique_ancestors.C` walks each positive-energy proton-beam-gas contribution through the first parent repeatedly until it reaches a root ancestor. It hashes that ancestor's PDG, status, four-momentum, energy, and vertex, then prints accepted events, beam-gas contributions, and unique ancestor hashes.

### Energy and MIP spectra

`lfhcal_energy_spectrum` sums active hit energy by event, decoded channel, and readout layer, then fills one channel-energy spectrum per layer.

`histogram_mip` under `lfhcal_analysis/muons` uses a single-muon file and:

```text
reject passive hits and invalid readout layers
require channel radius <= 1000 mm
sum hit contributions with individual energy >= 5e-4 GeV
accumulate energy by physical cell ID and by readout channel
fill one all-layer tile spectrum and seven channel spectra
```

`insert_histogram_mip` finds the first insert simulated-hit collection, accepts decoded layers 1 through 60, and fills the raw hit-energy spectrum from 0 to 0.03 GeV. The insert analyses use `4e-4 GeV` as their MIP constant.

### Geometry and bookkeeping

`channel_count` reads active `LFHCALHits`, stores each unique decoded channel in readout layers 0 through 6, and prints per-layer and total unique-channel counts.

`count_channels.C` recursively reads insert files, constructs unique 50 mm virtual channels separately for every physical layer and side, prints those counts, assigns them to the four side-0 or three side-1 angular regions, and prints channels per angular chip.

`plot_hits` stores one position per unique physical insert cell for layers 1 through 4 and sides 0 and 1. Its ROOT canvases overlay physical cells, the 50 mm virtual grid, angular chip boundaries, the side division, and the 163 mm inner-ring cut.

`bkg_composition` sums insert-hit contribution energy separately for the seven generator-status families and fills one hit-level energy spectrum per represented family.

### Plot comparison and table macros

`scripts/overlay.py` recursively finds matching canvas paths in two ROOT files, clones the first histogram or graph from each matching canvas, styles the pair, and writes comparison canvases while preserving directory paths and logarithmic-axis settings.

The LFHCAL `print_occupancy_tables.C` reads the `0.1`, `0.3`, and `0.5 MIP` graph points and writes channel p95/p99, chip p95/p99, and total detector data rate to `lfhcal.csv`. It derives displayed channel data rates from `hit_rate * 32 * 4`.

The insert `print_occupancy_tables.C` reads the `0.1`, `0.5`, and `1.5 MIP` graph points for all three layouts and both ring selections. It writes one CSV per configuration containing channel p95/p99, chip p95/p99, and left/right total data rates.

`overhead_estimate.C` is an alternate packing estimate at `0.5 * 4e-4 GeV`. It compares original cells and 50 mm tiles with both ring selections. Per event and side, it packs the active-channel count into groups of up to 36 and reports raw channel payload, packed-chip payload, and the conservative case that charges one framing term per channel. This macro uses count-based groups of 36, while the consolidated occupancy executable uses the geometric chip mappings described above.
