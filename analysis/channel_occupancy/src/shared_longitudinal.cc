#include "shared_longitudinal.h"

#include <TH1D.h>

#include <string>

namespace rates::channel_occupancy::longitudinal {
// ----------------------------------------------------------------------------------
// Helpers
// ----------------------------------------------------------------------------------
bool in_x_slice(const rates::LFHCALCellPosition& position) {
  return position.x_mm >= kXMinMM && position.x_mm < kXMaxMM;
}

// ----------------------------------------------------------------------------------
// Initialization
// ----------------------------------------------------------------------------------
void init_layers(std::vector<LayerAccum>& layers) {
  layers.assign(kNLayers, {});
  for (int layer = 0; layer < kNLayers; ++layer) {
    layers[layer].h_hits_evt = new TH1D(
        ("h_hits_evt_layer" + std::to_string(layer)).c_str(),
        ("LFHCAL layer " + std::to_string(layer) + ";hits/event;Events").c_str(),
        200, 0, 200);
  }
}

}  // namespace rates::channel_occupancy::longitudinal
