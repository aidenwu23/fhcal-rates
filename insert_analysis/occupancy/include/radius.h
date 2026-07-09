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
struct ThresholdAccum {
  std::unordered_map<rates::VirtualLFHCALChipID, double, rates::VirtualLFHCALChipIDHash> bits;
};

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void accumulate_event(std::array<ThresholdAccum, 3>& products,
                      const mip::EventEnergyMap& event_energy,
                      const rates::InsertToLFHCALMapper& mapper);

void write_output(TFile& output,
                  const std::array<ThresholdAccum, 3>& products,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::radius
