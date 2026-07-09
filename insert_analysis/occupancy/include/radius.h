#pragma once

#include "insert_to_lfhcal.h"
#include "mip.h"

#include <TFile.h>

#include <array>
#include <unordered_map>

namespace rates::insert_occupancy::radius {

// ----------------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------------
struct ThresholdSum {
  std::unordered_map<rates::VirtualLFHCALChipID, double, rates::VirtualLFHCALChipIDHash> payload_bits;
};

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<ThresholdSum, 3>& threshold_sums,
                      const mip::EventEnergyMap& channel_energy_sum,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TFile& output,
                  const std::array<ThresholdSum, 3>& threshold_sums,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::radius
