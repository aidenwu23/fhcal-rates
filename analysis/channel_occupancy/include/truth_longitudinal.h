#pragma once

#include <podio/Frame.h>

#include "shared_longitudinal.h"

#include <string>
#include <vector>

namespace rates::channel_occupancy::longitudinal {
// ----------------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------------
struct TruthGroup {
  std::string label;
  std::vector<LayerAccum> layers;
};

// ----------------------------------------------------------------------------------
// API
// ----------------------------------------------------------------------------------
void init_truth_groups(std::vector<TruthGroup>& groups);
bool process_truth_event(const podio::Frame& frame,
                         const rates::LFHCALCellIDDecoder& decoder,
                         const ThresholdsByLayer& thresholds_geV,
                         std::vector<TruthGroup>& groups);

}  // namespace rates::channel_occupancy::longitudinal
