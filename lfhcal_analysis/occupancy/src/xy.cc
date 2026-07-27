#include "occupancy/include/xy.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>

#include <array>
#include <cstdint>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace rates::lfhcal_analysis::occupancy::xy {
namespace {

constexpr std::array<const char*, kNOriginGroups> kOriginNames = {
    "DIS", "synrad", "eBrem", "eTouschek", "eCoulomb", "pBeamGas", "other", "all"};
constexpr double kYZSliceMinXMM = -400.0;
constexpr double kYZSliceMaxXMM = 300.0;

struct AxisEdges {
  std::vector<double> x;
  std::vector<double> y;
};

struct YBinStats {
  std::uint64_t passing_events = 0;
  std::size_t channels = 0;
};

std::vector<double> make_edges(const std::set<double>& coordinates) {
  std::vector<double> values(coordinates.begin(), coordinates.end());
  if (values.empty()) return {};
  if (values.size() == 1) return {values.front() - 0.5, values.front() + 0.5};

  std::vector<double> edges(values.size() + 1);
  edges.front() = values.front() - 0.5 * (values[1] - values[0]);
  for (std::size_t index = 1; index < values.size(); ++index) {
    edges[index] = 0.5 * (values[index - 1] + values[index]);
  }
  edges.back() = values.back() + 0.5 * (values.back() - values[values.size() - 2]);
  return edges;
}

AxisEdges spatial_edges(const LayerSum& layer_sum) {
  std::set<double> x_coordinates;
  std::set<double> y_coordinates;
  for (const auto& [channel, stats] : layer_sum.channels) {
    (void)channel;
    x_coordinates.insert(stats.x_mm);
    y_coordinates.insert(stats.y_mm);
  }
  return {make_edges(x_coordinates), make_edges(y_coordinates)};
}

void draw_and_write_spatial(TDirectory* directory, TH1* histogram, const char* canvas_name, bool two_dimensional) {
  directory->cd();
  TCanvas canvas(canvas_name, histogram->GetTitle(), 1000, 800);
  if (two_dimensional) canvas.SetLogz();
  histogram->SetStats(false);
  histogram->Draw(two_dimensional ? "colz" : "hist");
  rates::pad_axes(*histogram);
  canvas.Write();
}

void write_xy(TDirectory* parent, const Outputs& outputs, std::uint64_t n_events) {
  // Write aligned XY occupancy and rate maps for every origin and layer group.
  auto* xy_directory = parent->mkdir("xy");
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;

  // Use inclusive geometry to keep every origin map on common axes.
  for (std::size_t origin = 0; origin < kNOriginGroups; ++origin) {
    auto* origin_directory = xy_directory->mkdir(kOriginNames[origin]);

    // Write one transverse map per physical layer and one summed-layer map.
    for (int layer = 0; layer <= kAllLayers; ++layer) {
      const AxisEdges edges = spatial_edges(outputs.groups[kAllOrigins][layer]);
      if (edges.x.empty() || edges.y.empty()) continue;

      auto* layer_directory = origin_directory->mkdir(layer_name(layer).c_str());
      auto* hist_directory = layer_directory->mkdir("hists");
      const std::string title = std::string(kOriginNames[origin]) + " " + layer_title(layer);
      TH2D occupancy("h_occupancy",
                     (title + ";x [mm];y [mm];mean fired channel-events/event").c_str(),
                     static_cast<int>(edges.x.size()) - 1,
                     edges.x.data(),
                     static_cast<int>(edges.y.size()) - 1,
                     edges.y.data());
      TH2D rate("h_rate",
                (title + ";x [mm];y [mm];rate [Hz/channel]").c_str(),
                static_cast<int>(edges.x.size()) - 1,
                edges.x.data(),
                static_cast<int>(edges.y.size()) - 1,
                edges.y.data());

      // Fill the selected origin's channel values on the shared geometry.
      for (const auto& [channel, stats] : outputs.groups[origin][layer].channels) {
        (void)channel;
        occupancy.Fill(stats.x_mm,
                       stats.y_mm,
                       static_cast<double>(stats.passing_events) / static_cast<double>(n_events));
        rate.Fill(stats.x_mm,
                  stats.y_mm,
                  static_cast<double>(stats.passing_events) / total_time_sec);
      }

      hist_directory->cd();
      occupancy.Write();
      rate.Write();
      draw_and_write_spatial(layer_directory, &occupancy, "c_occupancy", true);
      draw_and_write_spatial(layer_directory, &rate, "c_rate", true);
    }
  }
}

void write_yz(TDirectory* parent, const Outputs& outputs, std::uint64_t n_events) {
  // Build a common Y axis from inclusive channels inside the longitudinal slice.
  std::set<double> y_coordinates;
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    for (const auto& [channel, stats] : outputs.groups[kAllOrigins][layer].channels) {
      (void)channel;
      if (stats.x_mm < kYZSliceMinXMM || stats.x_mm > kYZSliceMaxXMM) continue;
      y_coordinates.insert(stats.y_mm);
    }
  }
  const std::vector<double> y_edges = make_edges(y_coordinates);
  if (y_edges.empty()) return;

  auto* yz_directory = parent->mkdir("yz");

  // Collapse channels within each layer and Y row for every origin group.
  for (std::size_t origin = 0; origin < kNOriginGroups; ++origin) {
    std::array<std::unordered_map<double, YBinStats>, kNReadoutLayers> layer_bins;
    std::unordered_map<double, YBinStats> y_bins;

    // Accumulate channel pass totals in the central X slice.
    for (int layer = 0; layer < kNReadoutLayers; ++layer) {
      for (const auto& [channel, stats] : outputs.groups[origin][layer].channels) {
        (void)channel;
        if (stats.x_mm < kYZSliceMinXMM || stats.x_mm > kYZSliceMaxXMM) continue;
        auto& layer_bin = layer_bins[layer][stats.y_mm];
        layer_bin.passing_events += stats.passing_events;
        ++layer_bin.channels;
        auto& y_bin = y_bins[stats.y_mm];
        y_bin.passing_events += stats.passing_events;
        ++y_bin.channels;
      }
    }

    auto* origin_directory = yz_directory->mkdir(kOriginNames[origin]);
    auto* hist_directory = origin_directory->mkdir("hists");
    const std::string title = std::string(kOriginNames[origin]) + " central X slice";
    TH2D occupancy("h_occupancy",
                   (title + ";readout layer;y [mm];mean fired channel-events/event/channel").c_str(),
                   kNReadoutLayers,
                   -0.5,
                   static_cast<double>(kNReadoutLayers) - 0.5,
                   static_cast<int>(y_edges.size()) - 1,
                   y_edges.data());
    TH2D rate("h_rate",
              (title + ";readout layer;y [mm];rate [Hz/channel]").c_str(),
              kNReadoutLayers,
              -0.5,
              static_cast<double>(kNReadoutLayers) - 0.5,
              static_cast<int>(y_edges.size()) - 1,
              y_edges.data());
    TH1D y_projection("h_occupancy_y",
                      (title + ";y [mm];mean fired channel-events/event/channel").c_str(),
                      static_cast<int>(y_edges.size()) - 1,
                      y_edges.data());

    // Fill channel-averaged YZ bins and the collapsed Y projection.
    for (int layer = 0; layer < kNReadoutLayers; ++layer) {
      for (const auto& [y_mm, stats] : layer_bins[layer]) {
        const double mean_events = static_cast<double>(stats.passing_events) /
            (static_cast<double>(stats.channels) * static_cast<double>(n_events));
        occupancy.Fill(layer, y_mm, mean_events);
        rate.Fill(layer, y_mm, mean_events / kEventWindowSec);
      }
    }
    for (const auto& [y_mm, stats] : y_bins) {
      y_projection.Fill(y_mm,
                        static_cast<double>(stats.passing_events) /
                            (static_cast<double>(stats.channels) * static_cast<double>(n_events)));
    }

    hist_directory->cd();
    occupancy.Write();
    rate.Write();
    y_projection.Write();
    draw_and_write_spatial(origin_directory, &occupancy, "c_occupancy", true);
    draw_and_write_spatial(origin_directory, &rate, "c_rate", true);
    draw_and_write_spatial(origin_directory, &y_projection, "c_occupancy_y", false);
  }
}

}  // namespace

void accumulate_event(std::array<LayerSum, kNReadoutLayers + 1>& layer_sums,
                      std::array<TH1D*, kNChannelTypes>& channel_type_event_hists,
                      Outputs& outputs,
                      const EventChannels& passing_channels) {
  // Count passing channels per channel type.
  std::array<int, kNChannelTypes> channels_per_type{};

  // Accumulate inclusive and origin-separated channel maps.
  for (const auto& event_channel : passing_channels) {
    const int layer = event_channel.channel.rlayerz;
    if (layer < 0 || layer >= kNReadoutLayers) continue;

    // Fill the fixed-threshold channel-rate accumulators.
    for (int group : {layer, kAllLayers}) {
      auto& stats = layer_sums[group].channels[event_channel.channel];
      stats.x_mm = event_channel.x_mm;
      stats.y_mm = event_channel.y_mm;
      ++stats.passing_events;
      stats.total_hits += static_cast<std::uint64_t>(event_channel.hit_count);
    }
    ++channels_per_type[channel_type(layer)];

    // Fill the inclusive spatial group and every contributing origin group.
    for (std::size_t origin = 0; origin < kNOriginGroups; ++origin) {
      if (origin != kAllOrigins && event_channel.energy_by_origin[origin] <= 0.0) continue;
      for (int group : {layer, kAllLayers}) {
        auto& stats = outputs.groups[origin][group].channels[event_channel.channel];
        stats.x_mm = event_channel.x_mm;
        stats.y_mm = event_channel.y_mm;
        ++stats.passing_events;
        stats.total_hits += static_cast<std::uint64_t>(event_channel.hit_count);
      }
    }
  }

  // Record one event-wise sum for each channel type.
  for (int type = 0; type < kNChannelTypes; ++type) {
    channel_type_event_hists[type]->Fill(channels_per_type[type]);
  }
}

void write_output(TDirectory* parent, const Outputs& outputs, std::uint64_t n_events) {
  // Write transverse and longitudinal spatial products.
  write_xy(parent, outputs, n_events);
  write_yz(parent, outputs, n_events);
}

}  // namespace rates::lfhcal_analysis::occupancy::xy
