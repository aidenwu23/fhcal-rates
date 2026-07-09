#pragma once

#include <TDirectory.h>
#include <TH1.h>
#include <TH1D.h>

#include "insert_to_lfhcal.h"

#include <array>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <vector>

namespace rates::insert_occupancy {

constexpr int kNLayers = 60;
constexpr int kNSegments = 7;
constexpr int kAllSegmentsIndex = kNSegments;
constexpr double kEventWindowSec = 2e-6;
constexpr double kVirtualCellSizeMM = 50.0;
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;
constexpr std::array<int, kNSegments> kSegmentFirstLayers = {1, 6, 11, 21, 31, 41, 51};
constexpr std::array<int, kNSegments> kSegmentLastLayers = {5, 10, 20, 30, 40, 50, 60};

struct ChannelStats {
  double x_mm = 0.0;
  double y_mm = 0.0;
  std::uint64_t total_hits = 0;
};

struct LayerAccum {
  std::unordered_map<rates::VirtualLFHCALChannelID, ChannelStats, rates::VirtualLFHCALChannelIDHash> channels;
  TH1D* h_hits_evt = nullptr;
};

struct AxisEdges2D {
  std::vector<double> x_edges;
  std::vector<double> y_edges;
};

int segment_index(int layer);
const char* segment_dir_name(int segment);
std::string segment_title(int segment);
void draw_and_write(TDirectory* canvas_dir, TH1* hist, const char* canvas_name, bool logz = false, bool logy = false);
std::vector<double> make_axis_edges(const std::set<double>& coords);
AxisEdges2D make_layer_axis_edges(const LayerAccum& layer_accum);
void init_layer_accumulations(std::vector<LayerAccum>& layers);
void fill_event_histograms(std::vector<LayerAccum>& layers,
                           const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& event_hits);

}  // namespace rates::insert_occupancy
