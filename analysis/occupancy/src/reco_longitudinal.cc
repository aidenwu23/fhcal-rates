#include "reco_longitudinal.h"

#include <edm4eic/CalorimeterHitCollection.h>

#include "utils.h"

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace br::occupancy::longitudinal {
namespace {

constexpr const char* kRecoHitCollection = "LFHCALRecHits";

struct EventChannel {
  double energy_gev = 0.0;
  double y_mm = 0.0;
};

}  // namespace

void init_reco_layers(std::vector<LayerAccum>& layers) {
  init_layers(layers);
}

bool process_reco_event(const podio::Frame& frame,
                        const br::LFHCALCellIDDecoder& decoder,
                        const ThresholdsByLayer& thresholds_geV,
                        std::vector<LayerAccum>& layers) {
  if (!br::has_collection(frame, kRecoHitCollection)) return false;

  // Per-event channel counters.
  std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>> event_counts(kNLayers);
  // event_channels[layer][channel]: For this event, summed channel signal before applying threshold.
  std::vector<std::unordered_map<br::LFHCALChannelID, EventChannel, br::LFHCALChannelIDHash>> event_channels(kNLayers);
  std::vector<int> layer_totals(kNLayers, 0);

  // Grab and loop over hits.
  const auto& hits = frame.get<edm4eic::CalorimeterHitCollection>(kRecoHitCollection);
  for (const auto& hit : hits) {
    // Get layer.
    const int layer = hit.getLayer();
    if (layer < 0 || layer >= kNLayers) continue;

    // Decode cell ID and apply the x slice.
    const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
    const auto position = decoder.position(cell_id);
    if (!in_x_slice(position)) continue;

    // All channels with the same y value are stacked on top of each other.
    const auto channel_id = decoder.channel(cell_id);
    auto& event_channel = event_channels[layer][channel_id];
    event_channel.energy_gev += hit.getEnergy();
    event_channel.y_mm = position.y_mm;
  }

  for (int layer = 0; layer < kNLayers; ++layer) {
    for (const auto& [channel_id, event_channel] : event_channels[layer]) {
      // Apply summed channel threshold.
      if (event_channel.energy_gev <= thresholds_geV[layer]) continue;

      auto& stats = layers[layer].channels[channel_id];
      stats.y_mm = event_channel.y_mm;

      // Count this channel once for the current event after applying the summed channel threshold.
      ++event_counts[layer][channel_id];
      ++layer_totals[layer];
    }
  }

  // End of this event:
  // move the temporary per-event counts into the long-lived channel statistics.
  for (int layer = 0; layer < kNLayers; ++layer) {
    layers[layer].h_hits_evt->Fill(layer_totals[layer]);
    for (const auto& [channel_id, count] : event_counts[layer]) {
      auto& stats = layers[layer].channels[channel_id];

      // Add this event's contribution to the all-events total.
      stats.total_hits += static_cast<std::uint64_t>(count);

      // Keep the largest number of hits this channel ever saw in any single event.
      stats.max_hits_event = std::max(stats.max_hits_event, count);
    }
  }

  return true;
}

}  // namespace br::occupancy::longitudinal
