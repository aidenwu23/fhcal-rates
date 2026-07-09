#pragma once

#include "insert_to_lfhcal.h"

#include <TFile.h>

#include <array>
#include <cstdint>
#include <unordered_map>

namespace rates::insert_occupancy::mip {

// ----------------------------------------------------------------------------------
// Structs and aliases
// ----------------------------------------------------------------------------------
struct ChannelThresholdSum {
  std::unordered_map<rates::VirtualLFHCALChannelID, std::uint64_t, rates::VirtualLFHCALChannelIDHash> pass_counts;
};

struct ChipThresholdSum {
  std::unordered_map<rates::VirtualLFHCALChipID, std::uint64_t, rates::VirtualLFHCALChipIDHash> pass_counts;
};

struct DataThresholdSum {
  std::unordered_map<rates::VirtualLFHCALChipID, double, rates::VirtualLFHCALChipIDHash> payload_bits;
};

using EventEnergyMap = std::unordered_map<rates::VirtualLFHCALChannelID, double, rates::VirtualLFHCALChannelIDHash>;

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<ChannelThresholdSum, 16>& channel_threshold_sums,
                      std::array<ChipThresholdSum, 16>& chip_threshold_sums,
                      std::array<DataThresholdSum, 16>& data_threshold_sums,
                      const EventEnergyMap& channel_energies,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TFile& output,
                  const std::array<ChannelThresholdSum, 16>& channel_threshold_sums,
                  const std::array<ChipThresholdSum, 16>& chip_threshold_sums,
                  const std::array<DataThresholdSum, 16>& data_threshold_sums,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::mip
