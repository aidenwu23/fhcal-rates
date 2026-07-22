#pragma once

#include "occupancy/include/shared.h"

#include <array>

namespace rates::lfhcal_analysis::occupancy::xy {

void accumulate_event(std::array<LayerSum, kNReadoutLayers + 1>& layer_sums,
                      const EventChannels& passing_channels);

}  // namespace rates::lfhcal_analysis::occupancy::xy
