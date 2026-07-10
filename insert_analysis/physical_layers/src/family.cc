#include "physical_layers/include/family.h"

#include "physical_layers/include/channel.h"
#include "physical_layers/include/xy.h"

#include <unordered_map>

namespace rates::insert_analysis::physical_layers {

void init_outputs(Outputs& outputs) {
  // Allocate the per-layer and inclusive histograms for this output variant.
  init_segment_sums(outputs.segment_sums, outputs.mode);
}

void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper) {
  // Keep hit multiplicity and summed energy separately until thresholds are applied.
  std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash> channel_hit_counts;
  mip::EventEnergyMap channel_energies;

  // Map each physical hit into the virtual channel used by this layer layout.
  for (const auto& hit : event_hits) {
    auto channel = mapper.channel(hit.layer, hit.x_mm, hit.y_mm);
    channel.layer = mapped_layer(outputs.mode, hit.layer);
    if (channel.layer < 0) continue;
    if (!keep_channel(outputs.mode, channel)) continue;

    ++channel_hit_counts[channel]; // Keep multiplicity for occupancy products.
    channel_energies[channel] += hit.energy_gev; // Sum before applying the MIP threshold.
  }

  // Build the fixed 0.5-MIP channel view used by occupancy products.
  std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash> thresholded_channel_hit_counts;
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
      n_events);
}

}  // namespace rates::insert_analysis::physical_layers
