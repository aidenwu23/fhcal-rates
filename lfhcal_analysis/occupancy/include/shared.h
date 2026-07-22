#pragma once

#include <TDirectory.h>
#include <TH1.h>
#include <TH1D.h>

#include "decode_cell_id.h"
#include "event_channel.h"
#include "utils.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace rates::lfhcal_analysis::occupancy {

constexpr int kNReadoutLayers = 7;
constexpr int kAllLayers = kNReadoutLayers;
constexpr double kEventWindowSec = 2e-6;
constexpr double kNominalThresholdMIP = 0.1;
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;
constexpr std::array<double, 16> kMIPCoefficients = {
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
    0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

struct ChannelStats {
  // Full-sample position and nominal-threshold totals for one channel.
  double x_mm = 0.0;
  double y_mm = 0.0;
  std::uint64_t passing_events = 0;
  std::uint64_t total_hits = 0;
};

struct LayerSum {
  // Per-channel totals and event multiplicity for one layer group.
  std::unordered_map<rates::LFHCALChannelID, ChannelStats, rates::LFHCALChannelIDHash> channels;
  TH1D* h_channels_event = nullptr;
};

using EventChannels = std::vector<rates::lfhcal_analysis::EventChannel>;

std::string layer_name(int layer);
std::string layer_title(int layer);
std::size_t percentile_index(std::size_t n_values, double percentile);
void draw_and_write(TDirectory* directory, TH1* histogram, const char* canvas_name, bool logx = false, bool logy = false);

}  // namespace rates::lfhcal_analysis::occupancy
