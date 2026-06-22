#pragma once

#include <podio/Frame.h>

#include "shared.h"

#include <vector>

namespace br::occupancy {

void init_reco_layers(std::vector<LayerAccum>& layers);
bool process_reco_event(const podio::Frame& frame,
                        const br::LFHCALCellIDDecoder& decoder,
                        const ThresholdsByLayer& thresholds_geV,
                        std::vector<LayerAccum>& layers);

}  // namespace br::occupancy
