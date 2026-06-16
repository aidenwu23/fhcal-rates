#include "shared_longitudinal.h"

#include <TH1D.h>

#include <string>

namespace br::occupancy::longitudinal {

bool in_x_slice(const br::LFHCALCellPosition& position) {
  return position.x_mm >= kXMinMM && position.x_mm < kXMaxMM;
}

void init_layers(std::vector<LayerAccum>& layers) {
  layers.assign(kNLayers, {});
  for (int layer = 0; layer < kNLayers; ++layer) {
    layers[layer].h_hits_evt = new TH1D(
        ("h_hits_evt_layer" + std::to_string(layer)).c_str(),
        ("LFHCAL layer " + std::to_string(layer) + ";hits/event;Events").c_str(),
        200, 0, 200);
  }
}

}  // namespace br::occupancy::longitudinal
