# Insert analysis guide

`occupancy` reads simulated insert calorimeter hits and writes occupancy and
rate products for two longitudinal layouts.  Every accepted hit is first mapped onto
a 50 mm by 50 mm virtual LFHCAL channel.  The same event then feeds every requested
layout and inner-ring treatment.

## Main flow

1. `occupancy.cc` expands `-i` into ROOT files, finds the `HcalEndcapPInsert` hit
   collection, decodes each hit, and creates one `EventHit` per accepted hit.
2. It passes the same `EventHit` vector to four output accumulators:
   - physical layers with the inner ring;
   - physical layers with the inner ring removed;
   - LFHCAL-like longitudinal segments with the inner ring;
   - LFHCAL-like longitudinal segments with the inner ring removed.
3. Each accumulator maps hits to virtual channels, sums energy per channel within
   the event, applies the requested threshold, and adds the event to full-sample
   counters.
4. After every file is processed, each accumulator writes its products under its
   own directory in the output ROOT file.

`n_events` counts accepted events and converts accumulated counts or payload bits
into rates with the 2 microsecond event window.

## Shared input and mapping

| File | Role |
| --- | --- |
| `include/event_hit.h` | Small transport struct containing decoded layer, side, hit position in mm, and hit energy in GeV. |
| `scripts/include/decode_cell_id.h` and `scripts/src/decode_cell_id.cc` | Decode the insert cell ID into layer and side. |
| `scripts/include/insert_to_lfhcal.h` and `scripts/src/insert_to_lfhcal.cc` | Map `(layer, x, y)` to a 50 mm virtual channel, then group two by two channels into a virtual chip. |
| `scripts/include/utils.h` and `scripts/src/utils.cc` | Find ROOT files, check collections, and display file progress. |

`EventHit` remains independent of either output layout.  This keeps decoding in
`occupancy.cc` and leaves the layout-specific grouping in the two families below.

## Inner-ring variants

Each family has `with_inner_ring` and `no_inner_ring` wrapper files.  Their only job
is to initialize the output directory name and `OccupancyMode::drop_inner_ring`.
The generic family code receives that mode for every event.

`shared.cc::keep_channel` applies the removal.  The virtual-channel center is
measured from the global empty-hole center at `(-172 mm, 0 mm)`, combining the
insert placement and the hole offset inside the insert.  The removed disk has the
physical 146.1 mm beampipe-hole radius.  Since the filter runs before any occupancy,
chip, MIP, radius, or side product receives a channel, every statistic in a
`no_inner_ring` directory uses the same mask.

## Physical-layer family

Files under `physical_layers/` preserve all 60 physical layers.  A physical layer
maps to the zero-based channel layer index `physical_layer - 1`.

| File | Role |
| --- | --- |
| `include/family.h`, `src/family.cc` | Own `Outputs`, initialize its accumulators, process one event, and write one layout variant. |
| `include/shared.h`, `src/shared.cc` | Define common constants, channel statistics, output mode, layer names, inner-ring filtering, percentile indexing, and ROOT drawing helpers. |
| `include/xy.h`, `src/xy.cc` | Accumulate thresholded channel counts by layer and in the inclusive summed-layer entry. |
| `include/channel.h`, `src/channel.cc` | Write per-channel hit-rate and data-rate histograms for every layer and the layer sum. |
| `include/mip.h`, `src/mip.cc` | Scan 0.0 through 1.5 MIP thresholds and write channel hit-rate and data-rate percentile curves. |
| `include/with_inner_ring.h`, `src/with_inner_ring.cc` | Configure the `physical_layers/with_inner_ring` output variant. |
| `include/no_inner_ring.h`, `src/no_inner_ring.cc` | Configure the `physical_layers/no_inner_ring` output variant. |

For one event, `family.cc` keeps both hit multiplicity and summed channel energy.
The fixed 0.5-MIP selection feeds `xy.cc` and `channel.cc`; the unthresholded energy
map feeds the complete MIP scan in `mip.cc`.

## LFHCAL-segment family

Files under `lfhcal_segments/` combine physical layers into seven LFHCAL-like
longitudinal segments: `1-5`, `6-10`, `11-20`, `21-30`, `31-40`, `41-50`, and
`51-60`.  Every channel in a segment receives that segment's zero-based index as
its virtual-channel layer ID, so hits from all layers in that segment accumulate
together.

| File | Role |
| --- | --- |
| `include/family.h`, `src/family.cc` | Own all segment-level products and route one event into every product family. |
| `include/shared.h`, `src/shared.cc` | Define segment boundaries, map physical layers to segments, and provide shared filtering and ROOT helpers. |
| `include/xy.h`, `src/xy.cc` | Accumulate thresholded channel occupancy for each segment and the segment sum. |
| `include/channel.h`, `src/channel.cc` | Write virtual-channel hit-rate and data-rate histograms by segment. |
| `include/chip.h`, `src/chip.cc` | Combine thresholded virtual channels into virtual chips and write chip hit-rate and data-rate histograms. |
| `include/mip.h`, `src/mip.cc` | Produce MIP scans for both virtual channels and virtual chips. |
| `include/radius.h`, `src/radius.cc` | Bin thresholded virtual-chip data rates by radius and write p95 and p99 radius scans. |
| `include/side.h`, `src/side.cc` | Sum thresholded virtual-chip payload separately for the two insert sides and write the side data-rate scan. |
| `include/with_inner_ring.h`, `src/with_inner_ring.cc` | Configure the `lfhcal_segments/with_inner_ring` output variant. |
| `include/no_inner_ring.h`, `src/no_inner_ring.cc` | Configure the `lfhcal_segments/no_inner_ring` output variant. |

The segment family keeps three event-level views: channel hit multiplicity, summed
channel energy, and summed channel energy separated by insert side.  These supply
the occupancy, chip, MIP, radius, and side calculations without rereading hits.

## Other insert executables

| File | Executable | Purpose |
| --- | --- | --- |
| `histogram_mip.cc` | `insert_histogram_mip` | Fill one spectrum of accepted insert-hit energy, used to inspect the MIP response. |
| `plot_hits.cc` | `plot_hits` | Draw unique physical-cell positions for the first four insert layers, separated by side. |
| `bkg_composition.cc` | `bkg_composition` | Sum hit-contribution energy by generator-status origin and draw one spectrum per origin. |
| `CMakeLists.txt` | — | Builds these executables and `occupancy`, linking ROOT, podio, and EDM4hep. |

## Output tree

The occupancy ROOT file has two top-level layout directories:

```text
physical_layers/
  with_inner_ring/
  no_inner_ring/
lfhcal_segments/
  with_inner_ring/
  no_inner_ring/
```

Physical-layer variants contain channel occupancy and channel MIP products.
Segment variants contain those products plus virtual-chip, radius, and side products.
