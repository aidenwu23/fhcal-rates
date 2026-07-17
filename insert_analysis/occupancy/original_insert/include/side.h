#pragma once

#include "insert_to_lfhcal.h"
#include "original_insert/include/mip.h"

#include <TDirectory.h>

#include <array>

namespace rates::insert_analysis::original_insert::side {

struct ThresholdSum {
  double total_payload_bits = 0.0;
};

void accumulate_event(std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                      const OccupancyMode& mode,
                      const std::array<mip::EventEnergyMap, 2>& side_channel_energies,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TDirectory* parent,
                  const std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::original_insert::side
