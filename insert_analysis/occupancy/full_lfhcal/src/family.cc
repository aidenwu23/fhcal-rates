#include "full_lfhcal/include/family.h"

#include "full_lfhcal/include/channel.h"
#include "full_lfhcal/include/xy.h"

#include <unordered_map>

namespace rates::insert_analysis::full_lfhcal {

void init_outputs(Outputs& outputs) {
  // Allocate products shared by all events for this output variant.
  init_segment_sums(outputs.segment_sums, outputs.mode);
  chip::init_chip_sum(outputs.chip_sum, outputs.mode);
}

void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper) {
  // Keep hit multiplicity, total energy, and per-side energy until thresholds are applied.
  std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash> channel_hit_counts;
  mip::EventEnergyMap channel_energies;
  std::array<mip::EventEnergyMap, 2> side_channel_energies;

  // Map each physical hit into the virtual channel used by this segment layout.
  for (const auto& hit : event_hits) {
    auto channel = mapper.channel(hit.layer, hit.x_mm, hit.y_mm);
    channel.layer = mapped_layer(outputs.mode, hit.layer);
    if (channel.layer < 0) continue;

    // Keep or discard the channel depending on if we want to keep the inner ring.
    if (!keep_channel(outputs.mode, channel)) continue;

    ++channel_hit_counts[channel]; // Keep multiplicity for occupancy products.
    channel_energies[channel] += hit.energy_gev; // Sum energy before thresholding.
    side_channel_energies[hit.side][channel] += hit.energy_gev; // Keep each insert side separate.
  }

  // Build the fixed 0.5-MIP channel view used by occupancy products.
  std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash> thresholded_channel_hit_counts;
  for (const auto& [channel, energy] : channel_energies) {
    if (energy <= 0.5 * channel_mip_energy_gev(outputs.mode, channel.layer)) continue;
    auto it = channel_hit_counts.find(channel);
    if (it == channel_hit_counts.end()) continue;
    thresholded_channel_hit_counts[channel] = it->second;
  }

  // Send the thresholded view to occupancy products and the energy maps to scans.
  xy::accumulate_event(outputs.segment_sums, outputs.mode, thresholded_channel_hit_counts);
  chip::accumulate_event(outputs.chip_sum, thresholded_channel_hit_counts, mapper);
  mip::accumulate_event(
      outputs.channel_threshold_sums,
      outputs.channel_data_threshold_sums,
      outputs.chip_threshold_sums,
      outputs.chip_data_threshold_sums,
      outputs.mode,
      channel_energies,
      mapper);
  radius::accumulate_event(outputs.radius_threshold_sums, outputs.mode, channel_energies, mapper);
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
  radius::write_output(variant_dir, outputs.radius_threshold_sums, n_events);
  side::write_output(variant_dir, outputs.side_threshold_sums, outputs.mode, n_events);
}

}  // namespace rates::insert_analysis::full_lfhcal
