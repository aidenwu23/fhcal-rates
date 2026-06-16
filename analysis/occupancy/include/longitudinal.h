#pragma once

#include <podio/Frame.h>

#include "decode_cell_id.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

class TH1D;

namespace br::occupancy::longitudinal {

constexpr int kNLayers = 7;
constexpr double kEventWindowSec = 2e-6;
constexpr double kDisplayPaddingMM = 25.0;
constexpr double kXMinMM = -400.0;
constexpr double kXMaxMM = 300.0;

struct ChannelStats {
  double y_mm = 0.0;
  std::uint64_t total_hits = 0;
  int max_hits_event = 0;
};

struct LayerAccum {
  std::unordered_map<br::LFHCALChannelID, ChannelStats, br::LFHCALChannelIDHash> channels;
  TH1D* h_hits_evt = nullptr;
};

struct TruthGroup {
  std::string label;
  std::vector<LayerAccum> layers;
};

bool in_x_slice(const br::LFHCALCellPosition& position);
void init_layers(std::vector<LayerAccum>& layers);

void init_reco_layers(std::vector<LayerAccum>& layers);
bool process_reco_event(const podio::Frame& frame,
                        const br::LFHCALCellIDDecoder& decoder,
                        double threshold_geV,
                        std::vector<LayerAccum>& layers);

void init_truth_groups(std::vector<TruthGroup>& groups);
bool process_truth_event(const podio::Frame& frame,
                         const br::LFHCALCellIDDecoder& decoder,
                         double threshold_geV,
                         std::vector<TruthGroup>& groups);

}  // namespace br::occupancy::longitudinal
