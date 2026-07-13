#include "original_insert/include/family.h"

#include "original_insert/include/channel.h"
#include "original_insert/include/xy.h"

#include <unordered_map>

namespace rates::insert_analysis::original_insert {

void init_outputs(Outputs& outputs) {
  // Allocate the per-layer and inclusive histograms for this output variant.
  init_segment_sums(outputs.segment_sums, outputs.mode);
}

void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits) {
  // Keep hit multiplicity and summed energy separately until thresholds are applied.
  std::unordered_map<OriginalChannelID, int, OriginalChannelIDHash> channel_hit_counts;
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
    channel_energies[channel] += hit.energy_gev; // Sum before applying the MIP threshold.
    side_channel_energies[hit.side][channel] += hit.energy_gev;
  }

  // Build the fixed 0.5-MIP channel view used by occupancy products.
  std::unordered_map<OriginalChannelID, int, OriginalChannelIDHash> thresholded_channel_hit_counts;
  for (const auto& [channel, energy] : channel_energies) {
    if (energy <= 0.5 * channel_mip_energy_gev(outputs.mode, channel.layer)) continue;
    auto it = channel_hit_counts.find(channel);
    if (it == channel_hit_counts.end()) continue;
    thresholded_channel_hit_counts[channel] = it->second;
  }

  // Send the thresholded view to occupancy products and the full energy map to the scan.
  xy::accumulate_event(outputs.segment_sums, outputs.mode, thresholded_channel_hit_counts);
  mip::accumulate_event(
      outputs.channel_threshold_sums,
      outputs.channel_data_threshold_sums,
      outputs.mode,
      channel_energies);
  side::accumulate_event(outputs.side_threshold_sums, outputs.mode, side_channel_energies);
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
  mip::write_output(
      variant_dir,
      outputs.channel_threshold_sums,
      outputs.channel_data_threshold_sums,
      outputs.mode,
      n_events);
  side::write_output(variant_dir, outputs.side_threshold_sums, n_events);
}

}  // namespace rates::insert_analysis::original_insert
