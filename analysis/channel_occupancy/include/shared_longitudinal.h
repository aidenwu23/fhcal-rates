#pragma once

#include <TH1D.h>

#include "decode_cell_id.h"

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace rates::channel_occupancy::longitudinal {
// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
constexpr int kNLayers = 7;
constexpr double kEventWindowSec = 2e-6;
constexpr double kDisplayPaddingMM = 25.0;
constexpr double kXMinMM = -400.0;
constexpr double kXMaxMM = 300.0;
using ThresholdsByLayer = std::array<double, kNLayers>;

struct ChannelStats {
  double y_mm = 0.0;
  std::uint64_t total_hits = 0;
  int max_hits_event = 0;
};

struct LayerAccum {
  std::unordered_map<rates::LFHCALChannelID, ChannelStats, rates::LFHCALChannelIDHash> channels;
  TH1D* h_hits_evt = nullptr;
};

// ----------------------------------------------------------------------------------
// Helpers
// ----------------------------------------------------------------------------------
bool in_x_slice(const rates::LFHCALCellPosition& position);
void init_layers(std::vector<LayerAccum>& layers);

}  // namespace rates::channel_occupancy::longitudinal
