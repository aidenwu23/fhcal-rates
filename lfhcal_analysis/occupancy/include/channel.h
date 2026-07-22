#pragma once

#include "occupancy/include/shared.h"

#include <TDirectory.h>

#include <array>
#include <cstdint>

namespace rates::lfhcal_analysis::occupancy::channel {

void init(std::array<LayerSum, kNReadoutLayers + 1>& layer_sums);
void write_output(TDirectory* parent,
                  const std::array<LayerSum, kNReadoutLayers + 1>& layer_sums,
                  std::uint64_t n_events);

}  // namespace rates::lfhcal_analysis::occupancy::channel
