#include "occupancy/include/xy.h"

namespace rates::lfhcal_analysis::occupancy::xy {

void accumulate_event(std::array<LayerSum, kNReadoutLayers + 1>& layer_sums,
                      const EventChannels& passing_channels) {
  // Count passing channels per layer and across all layers.
  std::array<int, kNReadoutLayers + 1> channels_per_layer{};

  // Accumulate each passing channel into its layer and the inclusive layer sum.
  for (const auto& event_channel : passing_channels) {
    const int layer = event_channel.channel.rlayerz;
    if (layer < 0 || layer >= kNReadoutLayers) continue;

    // Fill the physical-layer and inclusive channel entries.
    for (int group : {layer, kAllLayers}) {
      auto& stats = layer_sums[group].channels[event_channel.channel];
      stats.x_mm = event_channel.x_mm;
      stats.y_mm = event_channel.y_mm;
      ++stats.passing_events;
      stats.total_hits += static_cast<std::uint64_t>(event_channel.hit_count);
      ++channels_per_layer[group];
    }
  }

  // Record this event's channel count in every layer group.
  for (int layer = 0; layer <= kAllLayers; ++layer) {
    layer_sums[layer].h_channels_event->Fill(channels_per_layer[layer]);
  }
}

}  // namespace rates::lfhcal_analysis::occupancy::xy
