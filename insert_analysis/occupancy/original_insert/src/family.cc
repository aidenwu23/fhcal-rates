#include "original_insert/include/family.h"

#include "original_insert/include/channel.h"
#include "original_insert/include/xy.h"

#include <unordered_map>

namespace rates::insert_analysis::original_insert {

void init_outputs(Outputs& outputs) {
  // Allocate the per-layer and inclusive histograms for this output variant.
  init_segment_sums(outputs.segment_sums, outputs.mode);
  chip::init_chip_sum(outputs.chip_sum, outputs.mode);
}

void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper) {
  // Keep hit multiplicity and summed energy separately until thresholds are applied.
  std::unordered_map<OriginalChannelID, int, OriginalChannelIDHash> channel_hit_counts;
  std::array<chip::ChannelHitCountMap, 2> side_channel_hit_counts;
  mip::EventEnergyMap channel_energies;
  std::array<mip::EventEnergyMap, 2> side_channel_energies;

  // Keep each original insert cell as its own channel in its physical layer.
  for (const auto& hit : event_hits) {
    const OriginalChannelID channel{
        hit.cell_id,
        mapped_layer(outputs.mode, hit.layer),
        hit.x_mm,
        hit.y_mm,
    };
    if (channel.layer < 0) continue;

    // Keep or discard the channel depending on if we want to keep the inner ring.
    if (!keep_channel(outputs.mode, channel)) continue;

    ++channel_hit_counts[channel]; // Keep multiplicity for occupancy products.
    ++side_channel_hit_counts[hit.side][channel];
    channel_energies[channel] += hit.energy_gev; // Sum before applying the MIP threshold.
    side_channel_energies[hit.side][channel] += hit.energy_gev;
  }

  // Build one thresholded channel map for general occupancy and two side-specific maps for chip assignment.
  std::unordered_map<OriginalChannelID, int, OriginalChannelIDHash> thresholded_channel_hit_counts;
  std::array<chip::ChannelHitCountMap, 2> thresholded_side_channel_hit_counts;

  // Loop thru all stored channels.
  for (const auto& [channel, energy] : channel_energies) {

    // Apply threshold.
    if (energy <= 0.5 * channel_mip_energy_gev(outputs.mode, channel.layer)) continue;

    // Recover this passing channel's hit multiplicity and add it to the general occupancy view.
    auto it = channel_hit_counts.find(channel);
    if (it == channel_hit_counts.end()) continue;
    thresholded_channel_hit_counts[channel] = it->second;

    for (int side = 0; side < 2; ++side) {
      const auto side_it = side_channel_hit_counts[side].find(channel);
      if (side_it != side_channel_hit_counts[side].end()) {
        // Copy the same passing channel and its side-specific hit multiplicity into the chip input map.
        thresholded_side_channel_hit_counts[side][channel] = side_it->second;
      }
    }
  }

  // Send the thresholded view to occupancy products and the full energy map to the scan.
  xy::accumulate_event(outputs.segment_sums, outputs.mode, thresholded_channel_hit_counts);
  chip::accumulate_event(outputs.chip_sum, thresholded_side_channel_hit_counts, mapper);
  mip::accumulate_event(
      outputs.channel_threshold_sums,
      outputs.channel_data_threshold_sums,
      outputs.chip_threshold_sums,
      outputs.chip_data_threshold_sums,
      outputs.mode,
      channel_energies,
      side_channel_energies,
      mapper);
  side::accumulate_event(outputs.side_threshold_sums, outputs.mode, side_channel_energies, mapper);
  fill_event_histograms(outputs.segment_sums, outputs.mode, thresholded_channel_hit_counts);
}

void write_output(TFile& output,
                  const Outputs& outputs,
                  std::uint64_t n_events) {
  // Reuse the layout directory while writing its separate inner-ring variants.
  auto* top_dir = output.GetDirectory(outputs.mode.top_dir_name);
  if (top_dir == nullptr) top_dir = output.mkdir(outputs.mode.top_dir_name);

  auto* variant_dir = top_dir->mkdir(outputs.mode.variant_dir_name);
  channel::write_output(variant_dir, outputs.mode, outputs.segment_sums, n_events);
  chip::write_output(variant_dir, outputs.mode, outputs.chip_sum, n_events);
  mip::write_output(
      variant_dir,
      outputs.channel_threshold_sums,
      outputs.channel_data_threshold_sums,
      outputs.chip_threshold_sums,
      outputs.chip_data_threshold_sums,
      outputs.mode,
      n_events);
  side::write_output(variant_dir, outputs.side_threshold_sums, n_events);
}

}  // namespace rates::insert_analysis::original_insert
