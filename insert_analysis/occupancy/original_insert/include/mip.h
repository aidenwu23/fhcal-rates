#pragma once

#include "original_insert/include/shared.h"

#include <TDirectory.h>

#include <array>
#include <unordered_map>

namespace rates::insert_analysis::original_insert::mip {

// ----------------------------------------------------------------------------------
// Structs and aliases
// ----------------------------------------------------------------------------------
struct ChannelThresholdSum {
  // pass_counts[channel]: events where the channel passed this threshold.
  std::unordered_map<OriginalChannelID, std::uint64_t, OriginalChannelIDHash> pass_counts;
};

struct ChannelDataThresholdSum {
  // payload_bits[channel]: data carried by the channel across the full sample.
  std::unordered_map<OriginalChannelID, double, OriginalChannelIDHash> payload_bits;
};

// EventEnergyMap[channel]: energy summed before thresholding one event.
using EventEnergyMap = std::unordered_map<OriginalChannelID, double, OriginalChannelIDHash>;

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<ChannelThresholdSum, 16>& channel_threshold_sums,
                      std::array<ChannelDataThresholdSum, 16>& channel_data_threshold_sums,
                      const OccupancyMode& mode,
                      const EventEnergyMap& channel_energies);

void write_output(TDirectory* parent,
                  const std::array<ChannelThresholdSum, 16>& channel_threshold_sums,
                  const std::array<ChannelDataThresholdSum, 16>& channel_data_threshold_sums,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::original_insert::mip
