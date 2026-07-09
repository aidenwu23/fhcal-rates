#pragma once

#include "shared.h"

#include <TFile.h>

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace rates::insert_occupancy::xy {

void accumulate_event(std::vector<LayerAccum>& layers,
                      const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& event_hits);

void write_output(TFile& output,
                  const std::vector<LayerAccum>& layers,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::xy
