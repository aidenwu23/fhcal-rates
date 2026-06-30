#pragma once

#include <podio/Frame.h>

#include "shared.h"
#include "decode_cell_id.h"

#include <string>
#include <vector>

namespace rates::channel_occupancy {
// ----------------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------------
struct TruthOccupancyGroup {
  std::string label;
  std::vector<LayerAccum> layers;
};

// ----------------------------------------------------------------------------------
// API
// ----------------------------------------------------------------------------------
void init_truth_groups(std::vector<TruthOccupancyGroup>& groups);
bool process_truth_event(const podio::Frame& frame,
                         const rates::LFHCALCellIDDecoder& decoder,
                         const ThresholdsByLayer& thresholds_geV,
                         std::vector<TruthOccupancyGroup>& groups);

}  // namespace rates::channel_occupancy
