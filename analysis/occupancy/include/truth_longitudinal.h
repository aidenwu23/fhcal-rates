#pragma once

#include <podio/Frame.h>

#include "shared_longitudinal.h"

#include <string>
#include <vector>

namespace br::occupancy::longitudinal {

struct TruthGroup {
  std::string label;
  std::vector<LayerAccum> layers;
};

void init_truth_groups(std::vector<TruthGroup>& groups);
bool process_truth_event(const podio::Frame& frame,
                         const br::LFHCALCellIDDecoder& decoder,
                         double threshold_geV,
                         std::vector<TruthGroup>& groups);

}  // namespace br::occupancy::longitudinal
