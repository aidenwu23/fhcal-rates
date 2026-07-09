#pragma once

#include "mip.h"
#include "insert_to_lfhcal.h"

#include <TFile.h>

#include <array>

namespace rates::insert_occupancy::side {

// ----------------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------------
struct ThresholdSum {
  double total_payload_bits = 0.0;
};

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                      const std::array<mip::EventEnergyMap, 2>& side_channel_energy_sum,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TFile& output,
                  const std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::side
