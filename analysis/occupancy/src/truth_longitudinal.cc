#include "longitudinal.h"

#include <TH1D.h>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "classify_hit.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace br::occupancy::longitudinal {
namespace {

constexpr const char* kTruthHitCollection = "LFHCALHits";

constexpr const char* kTruthLabels[] = {
    "DIS",
    "synrad",
    "eBrem",
    "eTouschek",
    "eCoulomb",
    "pBeamGas",
    "other",
    "all",
};

constexpr int kNTruthGroups = sizeof(kTruthLabels) / sizeof(kTruthLabels[0]);
constexpr int kAllTruthIndex = kNTruthGroups - 1;

}  // namespace

void init_truth_groups(std::vector<TruthGroup>& groups) {
  groups.clear();
  groups.reserve(kNTruthGroups);
  for (const char* label : kTruthLabels) {
    TruthGroup group;
    group.label = label;
    init_layers(group.layers);
    groups.push_back(std::move(group));
  }
}

bool process_truth_event(const podio::Frame& frame,
                         const br::LFHCALCellIDDecoder& decoder,
                         double threshold_geV,
                         std::vector<TruthGroup>& groups) {
  // Try to grab LFHCALHits.
  if (!br::has_collection(frame, kTruthHitCollection)) return false;

  // Per-event channel counters for every truth group and layer.
  std::vector<std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>>> event_counts(
      groups.size(), std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>>(kNLayers));
  std::vector<std::vector<int>> layer_totals(groups.size(), std::vector<int>(kNLayers, 0));

  // Grab and loop over hits.
  const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kTruthHitCollection);
  for (const auto& hit : hits) {
    // Apply threshold.
    if (hit.getEnergy() <= threshold_geV) continue;

    // Decode cell ID and apply the x slice.
    const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
    const auto position = decoder.position(cell_id);
    if (!in_x_slice(position)) continue;

    // Get readout layer and channel ID for this hit.
    const auto channel_id = decoder.channel(cell_id);
    const int layer = channel_id.rlayerz;
    if (layer < 0 || layer >= kNLayers) continue;

    std::array<double, kAllTruthIndex> energy_by_origin{};

    // Per hit, sum all contributions that fall into the same broad origin family.
    for (const auto& contribution : hit.getContributions()) {
      const double energy = contribution.getEnergy();
      if (energy <= 0.0) continue;

      const int origin_index = br::origin_index(contribution.getParticle().getGeneratorStatus());
      if (origin_index < 0 || origin_index >= kAllTruthIndex) continue;
      energy_by_origin[origin_index] += energy;
    }

    // Loop through all seen origins for this hit.
    for (int origin_index = 0; origin_index < kAllTruthIndex; ++origin_index) {
      // Skip if this origin contributes no energy.
      if (energy_by_origin[origin_index] <= 0.0) continue;

      // Otherwise fill the corresponding longitudinal bucket for this origin.
      auto& stats = groups[origin_index].layers[layer].channels[channel_id];
      stats.y_mm = position.y_mm;
      ++event_counts[origin_index][layer][channel_id];
      ++layer_totals[origin_index][layer];
    }

    // Also update the inclusive truth bucket once per hit.
    auto& all_stats = groups[kAllTruthIndex].layers[layer].channels[channel_id];
    all_stats.y_mm = position.y_mm;
    ++event_counts[kAllTruthIndex][layer][channel_id];
    ++layer_totals[kAllTruthIndex][layer];
  }

  // For each generatorStatus family...
  for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
    // For each readout layer...
    for (int layer = 0; layer < kNLayers; ++layer) {
      groups[group_index].layers[layer].h_hits_evt->Fill(layer_totals[group_index][layer]);

      // For each channel, accumulate stats.
      for (const auto& [channel_id, count] : event_counts[group_index][layer]) {
        auto& stats = groups[group_index].layers[layer].channels[channel_id];

        // Add this event's contribution to the all-events total.
        stats.total_hits += static_cast<std::uint64_t>(count);

        // Keep track of the single-event maximum for this channel.
        stats.max_hits_event = std::max(stats.max_hits_event, count);
      }
    }
  }

  return true;
}

}  // namespace br::occupancy::longitudinal
