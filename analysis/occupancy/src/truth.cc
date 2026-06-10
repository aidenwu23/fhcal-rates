#include "truth.h"

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "classify_hit.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <string>
#include <vector>

namespace br::occupancy {
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

void init_truth_groups(std::vector<TruthOccupancyGroup>& groups) {
  groups.clear();
  groups.reserve(kNTruthGroups);

  for (const char* label : kTruthLabels) {
    TruthOccupancyGroup group;
    group.label = label;
    init_reco_layers(group.layers);
    groups.push_back(std::move(group));
  }
}

bool process_truth_event(const podio::Frame& frame,
                         const br::LFHCALCellIDDecoder& decoder,
                         double threshold_geV,
                         std::vector<TruthOccupancyGroup>& groups) {
  if (!br::has_collection(frame, kTruthHitCollection)) return false;

  std::vector<std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>>> event_counts(
      groups.size(), std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>>(kNLayers + 1));
  std::vector<std::vector<int>> layer_totals(groups.size(), std::vector<int>(kNLayers + 1, 0));

  // Loop hits.
  const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kTruthHitCollection);
  for (const auto& hit : hits) {

    // Apply hit threshold.
    if (hit.getEnergy() <= threshold_geV) continue;

    // Decode cell_id.
    const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
    const auto channel_id = decoder.channel(cell_id);
    const int layer = channel_id.rlayerz;
    if (layer < 0 || layer >= kNLayers) continue;
    const auto cell_position = decoder.position(cell_id);

    std::array<double, kAllTruthIndex> energy_by_origin{};

    // Per hit, sum all contributions that fall into the same broad origin family.
    for (const auto& contribution : hit.getContributions()) {
      const double contribution_energy = contribution.getEnergy();
      if (contribution_energy <= 0.0) continue;

      const int origin_index = br::origin_index(contribution.getParticle().getGeneratorStatus());
      if (origin_index < 0 || origin_index >= kAllTruthIndex) continue;
      energy_by_origin[origin_index] += contribution_energy;
    }

    // Loop through all seen origins for this hit.
    for (int origin_index = 0; origin_index < kAllTruthIndex; ++origin_index) {

      // Skip if this origin contributes no energy.
      if (energy_by_origin[origin_index] <= 0.0) continue;

      // Otherwise fill the corresponding readout layer's plot for this origin once.
      auto fill_origin_bucket = [&](TruthOccupancyGroup& origin_bucket,
                                    int origin_bucket_index,
                                    int bucket_layer) {
        auto& stats = origin_bucket.layers[bucket_layer].channels[channel_id];
        stats.raw_cell_ids.insert(cell_id);
        stats.set_position(cell_position.x_mm, cell_position.y_mm);

        ++event_counts[origin_bucket_index][bucket_layer][channel_id];
        ++layer_totals[origin_bucket_index][bucket_layer];
      };

      fill_origin_bucket(groups[origin_index], origin_index, layer);
      fill_origin_bucket(groups[origin_index], origin_index, kAllLayersIndex);
    }

    // Also update the inclusive truth bucket once per hit.
    auto fill_all_truth_bucket = [&](TruthOccupancyGroup& all_bucket,
                                     int all_bucket_index,
                                     int bucket_layer) {
      auto& stats = all_bucket.layers[bucket_layer].channels[channel_id];
      stats.raw_cell_ids.insert(cell_id);
      stats.set_position(cell_position.x_mm, cell_position.y_mm);

      ++event_counts[all_bucket_index][bucket_layer][channel_id];
      ++layer_totals[all_bucket_index][bucket_layer];
    };
    fill_all_truth_bucket(groups[kAllTruthIndex], kAllTruthIndex, layer);
    fill_all_truth_bucket(groups[kAllTruthIndex], kAllTruthIndex, kAllLayersIndex);
  }

  // End of this event:
  // move the temporary per-event counts into the long-lived channel statistics.
  for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
    auto& group = groups[group_index];
    for (int layer = 0; layer <= kNLayers; ++layer) {
      group.layers[layer].h_hits_evt->Fill(layer_totals[group_index][layer]);
      for (const auto& [channel_id, count] : event_counts[group_index][layer]) {
        auto& stats = group.layers[layer].channels[channel_id];
        stats.total_hits += static_cast<std::uint64_t>(count);
        stats.max_hits_event = std::max(stats.max_hits_event, count);
      }
    }
  }

  return true;
}

}  // namespace br::occupancy
