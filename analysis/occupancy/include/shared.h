#pragma once

#include <TDirectory.h>
#include <TH1.h>
#include <TH1D.h>

#include "decode_cell_id.h"

#include <array>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace br::occupancy {

constexpr int kNLayers = 7;
constexpr int kAllLayersIndex = kNLayers;
constexpr double kEventWindowSec = 2e-6;
constexpr double kDisplayPaddingMM = 25.0;
using ThresholdsByLayer = std::array<double, kNLayers>;

struct ChannelStats {
  std::unordered_set<std::uint64_t> raw_cell_ids;

  double x_mm = 0.0;
  double y_mm = 0.0;

  std::uint64_t total_hits = 0;
  int max_hits_event = 0;

  double r() const;
};

struct LayerAccum {
  std::unordered_map<br::LFHCALChannelID, ChannelStats, br::LFHCALChannelIDHash> channels;
  TH1D* h_hits_evt = nullptr;
};

struct AxisEdges2D {
  std::vector<double> x_edges;
  std::vector<double> y_edges;
};

void draw_and_write(TDirectory* canvas_dir, TH1* hist, const char* canvas_name, bool logz = false, bool logy = false);
std::vector<double> make_axis_edges(const std::set<double>& coords);
AxisEdges2D make_layer_axis_edges(const LayerAccum& layer_accum);
void init_layer_accumulations(std::vector<LayerAccum>& layers);

}  // namespace br::occupancy
