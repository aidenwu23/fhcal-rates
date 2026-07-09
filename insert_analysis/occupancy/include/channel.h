#pragma once

#include "shared.h"

#include <TFile.h>

#include <cstdint>
#include <vector>

namespace rates::insert_occupancy::channel {

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void write_output(TFile& output,
                  const std::vector<LayerAccum>& layers,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::channel
