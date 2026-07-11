#pragma once

#include <TDirectory.h>
#include <TH1.h>
#include <TH1D.h>

#include "insert_to_lfhcal.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <vector>

namespace rates::insert_analysis::lfhcal_tiles {

constexpr int kNLayers = 60;
constexpr double kEventWindowSec = 2e-6;
constexpr double kVirtualCellSizeMM = 50.0;
constexpr double kInnerRingCenterXMM = -172.0;
constexpr double kInnerRingRadiusMM = 146.1;
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;
constexpr double kTileMipGeV = 4e-4;

struct OccupancyMode {
  // Directory names and channel selection for one output variant.
  const char* top_dir_name = "";
  const char* variant_dir_name = "";
  bool drop_inner_ring = false;
};

struct ChannelStats {
  // Channel center and hit total accumulated across all accepted events.
  double x_mm = 0.0;
  double y_mm = 0.0;
  std::uint64_t total_hits = 0;
  std::uint64_t passing_events = 0;
};

struct SegmentSum {
  // One layer, plus a separate inclusive layer-sum entry.
  std::unordered_map<rates::VirtualLFHCALChannelID, ChannelStats, rates::VirtualLFHCALChannelIDHash> channels;
  TH1D* h_hits_evt = nullptr;
};

struct AxisEdges2D {
  // Histogram edges built from the occupied virtual-channel centers.
  std::vector<double> x_edges;
  std::vector<double> y_edges;
};

int mapped_layer(const OccupancyMode& mode, int physical_layer);
int group_count(const OccupancyMode& mode);
int all_groups_index(const OccupancyMode& mode);
std::string group_dir_name(const OccupancyMode& mode, int group);
std::string group_title(const OccupancyMode& mode, int group);
double channel_mip_energy_gev(const OccupancyMode& mode, int mapped_layer);
bool keep_channel(const OccupancyMode& mode, const rates::VirtualLFHCALChannelID& channel);
std::size_t percentile_index(std::size_t n_values, double percentile);
void draw_and_write(TDirectory* canvas_dir, TH1* hist, const char* canvas_name, bool logz = false, bool logy = false);
std::vector<double> make_axis_edges(const std::set<double>& coords);
AxisEdges2D make_segment_axis_edges(const SegmentSum& segment_sum);
void init_segment_sums(std::vector<SegmentSum>& segment_sums, const OccupancyMode& mode);
void fill_event_histograms(std::vector<SegmentSum>& segment_sums,
                           const OccupancyMode& mode,
                           const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& channel_hit_counts);

}  // namespace rates::insert_analysis::lfhcal_tiles
