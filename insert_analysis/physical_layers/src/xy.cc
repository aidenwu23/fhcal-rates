#include "physical_layers/include/xy.h"

namespace rates::insert_analysis::physical_layers::xy {

// Called per event after the hit loop to accumulate x-y occupancy.
void accumulate_event(std::vector<SegmentSum>& segment_sums,
                      const OccupancyMode& mode,
                      const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& channel_hit_counts) {
  // Each event contributes to its own layer and to the inclusive summed-layer map.
  for (const auto& [channel, count] : channel_hit_counts) {
    if (channel.layer < 0 || channel.layer >= all_groups_index(mode)) continue;
    if (!keep_channel(mode, channel)) continue;

    // This lookup creates the channel entry on its first accepted event.
    auto& segment_sum = segment_sums[channel.layer];
    auto& stats = segment_sum.channels[channel];

    // Map the channel x and y index back onto the corresponding spatial location using cell size.
    stats.x_mm = (static_cast<double>(channel.ix) + 0.5) * kVirtualCellSizeMM;
    stats.y_mm = (static_cast<double>(channel.iy) + 0.5) * kVirtualCellSizeMM;
    stats.total_hits += count; // Add this event's accepted hit multiplicity.

    // Fill the inclusive layer sum with the same channel contribution.
    auto& all_segment_stats = segment_sums[all_groups_index(mode)].channels[channel];
    all_segment_stats.x_mm = stats.x_mm;
    all_segment_stats.y_mm = stats.y_mm;
    all_segment_stats.total_hits += count;
  }
}

}  // namespace rates::insert_analysis::physical_layers::xy
