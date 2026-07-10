#pragma once

#include "full_lfhcal/include/mip.h"
#include "insert_to_lfhcal.h"

#include <TDirectory.h>

#include <array>

namespace rates::insert_analysis::full_lfhcal::side {

// ----------------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------------
struct ThresholdSum {
  // Payload carried by this side across all accepted events.
  double total_payload_bits = 0.0;
};

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                      const OccupancyMode& mode,
                      const std::array<mip::EventEnergyMap, 2>& side_channel_energies,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TDirectory* parent,
                  const std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                  const OccupancyMode& mode,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::full_lfhcal::side
