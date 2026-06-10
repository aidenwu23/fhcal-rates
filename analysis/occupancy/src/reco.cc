#include "reco.h"

#include <edm4eic/CalorimeterHitCollection.h>

#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace br::occupancy {
namespace {

constexpr const char* kRecoHitCollection = "LFHCALRecHits";
}  // namespace

double ChannelStats::x() const { return n_pos > 0 ? x_sum / static_cast<double>(n_pos) : 0.0; }
double ChannelStats::y() const { return n_pos > 0 ? y_sum / static_cast<double>(n_pos) : 0.0; }
double ChannelStats::r() const { return std::hypot(x(), y()); }

void init_reco_layers(std::vector<LayerAccum>& layers) {
  layers.assign(kNLayers + 1, {});

  // Initialize layer accumulators and histograms.
  for (int layer = 0; layer <= kNLayers; ++layer) {
    layers[layer].h_hits_evt = new TH1D(
        "h_hits_evt",
        (layer == kAllLayersIndex ? std::string("LFHCAL summed layers;hits/event;Events")
          : std::string("LFHCAL layer ") + std::to_string(layer) + ";hits/event;Events").c_str(),
        200, 0, 200);
  }
}

bool process_reco_event(const podio::Frame& frame,
                        const br::LFHCALDecoder& decoder,
                        double threshold_geV,
                        std::vector<LayerAccum>& layers) {
  // Grab LFHCALRecHits when possible.
  if (!br::has_collection(frame, kRecoHitCollection)) return false;

  // Per-event channel counters.
  // event_counts[layer][channel]: For this event, the number of hits that landed in [layer]'s [channel].
  std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>> event_counts(kNLayers + 1);
  std::vector<int> layer_totals(kNLayers + 1, 0);

  // Grab and loop over hits.
  const auto& hits = frame.get<edm4eic::CalorimeterHitCollection>(kRecoHitCollection);
  for (const auto& hit : hits) {

    // Apply hit threshold.
    if (hit.getEnergy() <= threshold_geV) continue;

    // Get layer and cellID.
    const int layer = hit.getLayer();
    if (layer < 0 || layer >= kNLayers) continue;

    // Decode channel ID from cell ID.
    const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
    const auto channel_id = decoder.channel(cell_id);
    const auto pos = hit.getPosition();

    // Does one of two things:
    // 1. There are no preeixsting channels with this ID --> create new channel.
    // 2. There is a preexisting channel with this ID --> use that one.
    auto& stats = layers[layer].channels[channel_id];
    stats.raw_cell_ids.insert(cell_id);
    stats.x_sum += pos.x;
    stats.y_sum += pos.y;
    ++stats.n_pos;

    // Count one more hit in this channel for the current event.
    ++event_counts[layer][channel_id];
    ++layer_totals[layer];

    // Also update the "all layers together" bucket.
    auto& merged_stats = layers[kAllLayersIndex].channels[channel_id];
    merged_stats.raw_cell_ids.insert(cell_id);
    merged_stats.x_sum += pos.x;
    merged_stats.y_sum += pos.y;
    ++merged_stats.n_pos;
    ++event_counts[kAllLayersIndex][channel_id];
    ++layer_totals[kAllLayersIndex];
  }

  // End of this event:
  // move the temporary per-event counts into the long-lived channel statistics.
  for (int layer = 0; layer <= kNLayers; ++layer) {
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
}  // namespace br::occupancy
