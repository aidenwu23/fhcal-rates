#include "truth.h"

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "classify_hit.h"
#include "reco.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>

namespace rates::channel_occupancy {
namespace {
// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
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

struct EventChannel {
  double energy_gev = 0.0;
  std::array<double, kAllTruthIndex> energy_by_origin{};
  std::unordered_set<std::uint64_t> raw_cell_ids;
  double x_mm = 0.0;
  double y_mm = 0.0;
};

}  // namespace

// ----------------------------------------------------------------------------------
// Initialization
// ----------------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------------
// Event processing
// ----------------------------------------------------------------------------------
bool process_truth_event(const podio::Frame& frame,
                         const rates::LFHCALCellIDDecoder& decoder,
                         const ThresholdsByLayer& thresholds_geV,
                         std::vector<TruthOccupancyGroup>& groups) {
  if (!rates::has_collection(frame, kTruthHitCollection)) return false;

  std::vector<std::vector<std::unordered_map<rates::LFHCALChannelID, int, rates::LFHCALChannelIDHash>>> event_counts(
      groups.size(), std::vector<std::unordered_map<rates::LFHCALChannelID, int, rates::LFHCALChannelIDHash>>(kNLayers + 1));
  // event_channels[layer][channel]: For this event, summed channel signal before applying threshold.
  std::vector<std::unordered_map<rates::LFHCALChannelID, EventChannel, rates::LFHCALChannelIDHash>> event_channels(kNLayers);
  std::vector<std::vector<int>> layer_totals(groups.size(), std::vector<int>(kNLayers + 1, 0));

  // Loop hits.
  const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kTruthHitCollection);
  for (const auto& hit : hits) {

    // Decode cell_id.
    const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
    const auto channel_id = decoder.channel(cell_id);
    const int layer = channel_id.rlayerz;
    if (layer < 0 || layer >= kNLayers) continue;

    const auto cell_position = decoder.position(cell_id);

    // Does one of two things:
    // 1. There are no preexisting event channels with this ID --> create new event channel.
    // 2. There is a preexisting event channel with this ID --> accumulate into that one.
    auto& event_channel = event_channels[layer][channel_id];
    event_channel.energy_gev += hit.getEnergy(); // Sum energy here to apply a threshold later.
    event_channel.raw_cell_ids.insert(cell_id);
    event_channel.x_mm = cell_position.x_mm;
    event_channel.y_mm = cell_position.y_mm;

    // Per hit, sum all contributions that fall into the same broad origin family.
    for (const auto& contribution : hit.getContributions()) {
      const double contribution_energy = contribution.getEnergy();
      if (contribution_energy <= 0.0) continue;

      const int origin_index = rates::origin_index(contribution.getParticle().getGeneratorStatus());
      if (origin_index < 0 || origin_index >= kAllTruthIndex) continue;
      event_channel.energy_by_origin[origin_index] += contribution_energy;
    }
  }

  for (int layer = 0; layer < kNLayers; ++layer) {
    for (const auto& [channel_id, event_channel] : event_channels[layer]) {
      // Apply summed channel threshold.
      if (event_channel.energy_gev <= thresholds_geV[layer]) continue;

      // Loop through all seen origins for this channel.
      for (int origin_index = 0; origin_index < kAllTruthIndex; ++origin_index) {

        // Skip if this origin contributes no energy.
        if (event_channel.energy_by_origin[origin_index] <= 0.0) continue;

        // Otherwise fill the corresponding readout layer's plot for this origin.
        auto fill_origin_bucket = [&](TruthOccupancyGroup& origin_group,
                                      int group_index,
                                      int layer_index) {
          // Create/reference stats for this layer's channel.
          auto& stats = origin_group.layers[layer_index].channels[channel_id];
          stats.raw_cell_ids.insert(event_channel.raw_cell_ids.begin(), event_channel.raw_cell_ids.end());
          stats.x_mm = event_channel.x_mm;
          stats.y_mm = event_channel.y_mm;

          ++event_counts[group_index][layer_index][channel_id]; // Increment channel count for this origin.
          ++layer_totals[group_index][layer_index];
        };

        fill_origin_bucket(groups[origin_index], origin_index, layer);
        fill_origin_bucket(groups[origin_index], origin_index, kAllLayersIndex);
      }

      // Also update the inclusive truth bucket once per channel.
      auto fill_all_truth_bucket = [&](TruthOccupancyGroup& all_group,
                                       int all_group_index,
                                       int layer_index) {
        auto& stats = all_group.layers[layer_index].channels[channel_id];
        stats.raw_cell_ids.insert(event_channel.raw_cell_ids.begin(), event_channel.raw_cell_ids.end());
        stats.x_mm = event_channel.x_mm;
        stats.y_mm = event_channel.y_mm;

        ++event_counts[all_group_index][layer_index][channel_id];
        ++layer_totals[all_group_index][layer_index];
      };
      fill_all_truth_bucket(groups[kAllTruthIndex], kAllTruthIndex, layer);
      fill_all_truth_bucket(groups[kAllTruthIndex], kAllTruthIndex, kAllLayersIndex);
    }
  }

  // For each generatorStatus family...
  for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
    auto& group = groups[group_index];

    // For each readout layer...
    for (int layer = 0; layer <= kNLayers; ++layer) {
      group.layers[layer].h_hits_evt->Fill(layer_totals[group_index][layer]);

      // For each channel, acumulate stats.
      for (const auto& [channel_id, count] : event_counts[group_index][layer]) {
        auto& stats = group.layers[layer].channels[channel_id];
        stats.total_hits += static_cast<std::uint64_t>(count);  // Increment total hit count from this event.
        stats.max_hits_event = std::max(stats.max_hits_event, count); // Keep track of one with most hits.
      }
    }
  }

  return true;
}

}  // namespace rates::channel_occupancy
