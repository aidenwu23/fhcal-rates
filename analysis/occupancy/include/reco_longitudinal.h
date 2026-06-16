#pragma once

#include <podio/Frame.h>

#include "shared_longitudinal.h"

#include <vector>

namespace br::occupancy::longitudinal {

void init_reco_layers(std::vector<LayerAccum>& layers);
bool process_reco_event(const podio::Frame& frame,
                        const br::LFHCALCellIDDecoder& decoder,
                        double threshold_geV,
                        std::vector<LayerAccum>& layers);

}  // namespace br::occupancy::longitudinal
