#pragma once

#include "insert_to_lfhcal.h"
#include "full_lfhcal/include/shared.h"

#include <TDirectory.h>

#include <array>
#include <unordered_map>

namespace rates::insert_analysis::full_lfhcal::mip {

// ----------------------------------------------------------------------------------
// Structs and aliases
// ----------------------------------------------------------------------------------
struct ChannelThresholdSum {
  // pass_counts[channel]: events where the channel passed this threshold.
  std::unordered_map<rates::VirtualLFHCALChannelID, std::uint64_t, rates::VirtualLFHCALChannelIDHash> pass_counts;
};

struct ChannelDataThresholdSum {
  // payload_bits[channel]: data carried by the channel across the full sample.
  std::unordered_map<rates::VirtualLFHCALChannelID, double, rates::VirtualLFHCALChannelIDHash> payload_bits;
};

struct ChipThresholdSum {
  // pass_counts[chip]: events where at least one chip channel passed this threshold.
  std::unordered_map<rates::VirtualLFHCALChipID, std::uint64_t, rates::VirtualLFHCALChipIDHash> pass_counts;
};

struct ChipDataThresholdSum {
  // payload_bits[chip]: data carried by the chip across the full sample.
  std::unordered_map<rates::VirtualLFHCALChipID, double, rates::VirtualLFHCALChipIDHash> payload_bits;
};

// EventEnergyMap[channel]: energy summed before thresholding one event.
using EventEnergyMap = std::unordered_map<rates::VirtualLFHCALChannelID, double, rates::VirtualLFHCALChannelIDHash>;

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<ChannelThresholdSum, 16>& channel_threshold_sums,
                      std::array<ChannelDataThresholdSum, 16>& channel_data_threshold_sums,
                      std::array<ChipThresholdSum, 16>& chip_threshold_sums,
                      std::array<ChipDataThresholdSum, 16>& chip_data_threshold_sums,
                      const OccupancyMode& mode,
                      const EventEnergyMap& channel_energies,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TDirectory* parent,
                  const std::array<ChannelThresholdSum, 16>& channel_threshold_sums,
                  const std::array<ChannelDataThresholdSum, 16>& channel_data_threshold_sums,
                  const std::array<ChipThresholdSum, 16>& chip_threshold_sums,
                  const std::array<ChipDataThresholdSum, 16>& chip_data_threshold_sums,
                  const OccupancyMode& mode,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::full_lfhcal::mip
