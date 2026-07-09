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
struct ChannelThresholdAccum {
  std::unordered_map<rates::VirtualLFHCALChannelID, std::uint64_t, rates::VirtualLFHCALChannelIDHash> passes;
};

struct ChipThresholdAccum {
  std::unordered_map<rates::VirtualLFHCALChipID, std::uint64_t, rates::VirtualLFHCALChipIDHash> passes;
};

struct DataThresholdAccum {
  std::unordered_map<rates::VirtualLFHCALChipID, double, rates::VirtualLFHCALChipIDHash> bits;
};

using EventEnergyMap = std::unordered_map<rates::VirtualLFHCALChannelID, double, rates::VirtualLFHCALChannelIDHash>;

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<ChannelThresholdAccum, 16>& channel_products,
                      std::array<ChipThresholdAccum, 16>& chip_products,
                      std::array<DataThresholdAccum, 16>& data_products,
                      const EventEnergyMap& event_energy,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TFile& output,
                  const std::array<ChannelThresholdAccum, 16>& channel_products,
                  const std::array<ChipThresholdAccum, 16>& chip_products,
                  const std::array<DataThresholdAccum, 16>& data_products,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::mip
