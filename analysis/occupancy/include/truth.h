#pragma once

#include <podio/Frame.h>

#include "shared.h"
#include "decode_cell_id.h"

#include <string>
#include <vector>

namespace br::occupancy {

struct TruthOccupancyGroup {
  std::string label;
  std::vector<LayerAccum> layers;
};

void init_truth_groups(std::vector<TruthOccupancyGroup>& groups);
bool process_truth_event(const podio::Frame& frame,
                         const br::LFHCALCellIDDecoder& decoder,
                         const ThresholdsByLayer& thresholds_geV,
                         std::vector<TruthOccupancyGroup>& groups);

}  // namespace br::occupancy
