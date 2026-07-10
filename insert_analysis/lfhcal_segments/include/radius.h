#pragma once

#include "insert_to_lfhcal.h"
#include "lfhcal_segments/include/mip.h"

#include <TDirectory.h>

#include <array>
#include <unordered_map>

namespace rates::insert_analysis::lfhcal_segments::radius {

// ----------------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------------
struct ThresholdSum {
  // payload_bits[chip]: thresholded chip payload accumulated across the sample.
  std::unordered_map<rates::VirtualLFHCALChipID, double, rates::VirtualLFHCALChipIDHash> payload_bits;
};

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<ThresholdSum, 3>& threshold_sums,
                      const OccupancyMode& mode,
                      const mip::EventEnergyMap& channel_energies,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TDirectory* parent,
                  const std::array<ThresholdSum, 3>& threshold_sums,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::lfhcal_segments::radius
