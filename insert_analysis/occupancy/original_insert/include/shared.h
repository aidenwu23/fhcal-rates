#pragma once

#include <TDirectory.h>
#include <TH1.h>
#include <TH1D.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

namespace rates::insert_analysis::original_insert {

constexpr int kNLayers = 60;
constexpr double kEventWindowSec = 2e-6;
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
};

struct OriginalChannelID {
  std::uint64_t cell_id = 0;
  int layer = 0;
  double x_mm = 0.0;
  double y_mm = 0.0;

  bool operator==(const OriginalChannelID& other) const {
    return cell_id == other.cell_id;
  }
};

struct OriginalChannelIDHash {
  std::size_t operator()(const OriginalChannelID& channel) const {
    return std::hash<std::uint64_t>{}(channel.cell_id);
  }
};

struct SegmentSum {
  // One layer, plus a separate inclusive layer-sum entry.
  std::unordered_map<OriginalChannelID, ChannelStats, OriginalChannelIDHash> channels;
  TH1D* h_hits_evt = nullptr;
};

int mapped_layer(const OccupancyMode& mode, int physical_layer);
int group_count(const OccupancyMode& mode);
int all_groups_index(const OccupancyMode& mode);
std::string group_dir_name(const OccupancyMode& mode, int group);
std::string group_title(const OccupancyMode& mode, int group);
double channel_mip_energy_gev(const OccupancyMode& mode, int mapped_layer);
bool keep_channel(const OccupancyMode& mode, const OriginalChannelID& channel);
std::size_t percentile_index(std::size_t n_values, double percentile);
void draw_and_write(TDirectory* canvas_dir, TH1* hist, const char* canvas_name, bool logz = false, bool logy = false);
void init_segment_sums(std::vector<SegmentSum>& segment_sums, const OccupancyMode& mode);
void fill_event_histograms(std::vector<SegmentSum>& segment_sums,
                           const OccupancyMode& mode,
                           const std::unordered_map<OriginalChannelID, int, OriginalChannelIDHash>& channel_hit_counts);

}  // namespace rates::insert_analysis::original_insert
