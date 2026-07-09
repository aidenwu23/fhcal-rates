#pragma once

#include "mip.h"
#include "insert_to_lfhcal.h"

#include <TFile.h>

#include <array>

namespace rates::insert_occupancy::side {

// ----------------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------------
struct ThresholdAccum {
  double total_bits = 0.0;
};

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<std::array<ThresholdAccum, 16>, 2>& products,
                      const std::array<mip::EventEnergyMap, 2>& side_event_energy,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TFile& output,
                  const std::array<std::array<ThresholdAccum, 16>, 2>& products,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::side
