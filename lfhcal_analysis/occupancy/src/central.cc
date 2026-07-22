#include "occupancy/include/central.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1D.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace rates::lfhcal_analysis::occupancy::central {
namespace {

constexpr double kMaxRadiusMM = 1000.0;
constexpr std::array<double, 2> kPercentiles = {0.95, 0.99};
constexpr std::array<int, 2> kColors = {kBlue + 1, kRed + 1};
constexpr std::array<std::size_t, 3> kDistributionThresholds = {1, 5, 15};
constexpr std::array<const char*, 3> kSourceNames = {"all", "DIS", "pBeamGas"};

const ChannelPassMap& source_map(const Outputs& outputs, int source, std::size_t threshold) {
  return source == 0 ? outputs.all_passes[threshold] : outputs.origin_passes[source - 1][threshold];
}

std::array<std::array<double, kMIPCoefficients.size()>, 2> percentile_values(
    const Outputs& outputs,
    int source,
    int layer,
    double total_time_sec) {
  // Compute one p95 and p99 channel rate per threshold for this source and layer.
  std::array<std::array<double, kMIPCoefficients.size()>, 2> values{};

  // Loop over thresholds.
  for (std::size_t threshold = 0; threshold < kMIPCoefficients.size(); ++threshold) {
    std::vector<double> rates_hz;

    // Collect rates from channels in the requested layer.
    for (const auto& [channel, passes] : source_map(outputs, source, threshold)) {
      if (channel.rlayerz != layer) continue;
      rates_hz.push_back(static_cast<double>(passes) / total_time_sec);
    }
    if (rates_hz.empty()) continue;
    std::sort(rates_hz.begin(), rates_hz.end());
    for (std::size_t percentile = 0; percentile < kPercentiles.size(); ++percentile) {
      values[percentile][threshold] = rates_hz[percentile_index(rates_hz.size(), kPercentiles[percentile])];
    }
  }
  return values;
}

void write_graphs(TDirectory* directory,
                  const char* source,
                  int layer,
                  const std::array<std::array<double, kMIPCoefficients.size()>, 2>& values) {
  // Draw the p95 and p99 threshold scans for one source and layer.
  directory->cd();
  const std::string title = "LFHCAL " + std::string(source) + " channel rates, R < 1000 mm, layer " +
                            std::to_string(layer + 1) + ";MIP threshold;rate [Hz]";
  TCanvas canvas(("c_" + std::string(source) + "_layer" + std::to_string(layer + 1)).c_str(), title.c_str(), 1000, 800);
  canvas.SetGrid();
  auto* frame = canvas.DrawFrame(0.0, 0.0, 1.5, 2.0e5);
  frame->SetTitle(title.c_str());
  TLegend legend(0.7, 0.76, 0.88, 0.88);
  std::array<TGraph, 2> graphs = {
      TGraph(static_cast<int>(kMIPCoefficients.size()), kMIPCoefficients.data(), values[0].data()),
      TGraph(static_cast<int>(kMIPCoefficients.size()), kMIPCoefficients.data(), values[1].data())};

  // Style, label, and write both percentile graphs.
  for (std::size_t percentile = 0; percentile < graphs.size(); ++percentile) {
    const std::string label = "p" + std::to_string(static_cast<int>(100.0 * kPercentiles[percentile]));
    graphs[percentile].SetName(("g_" + std::string(source) + "_layer" + std::to_string(layer + 1) + "_" + label).c_str());
    graphs[percentile].SetLineColor(kColors[percentile]);
    graphs[percentile].SetMarkerColor(kColors[percentile]);
    graphs[percentile].SetMarkerStyle(20 + static_cast<int>(percentile));
    graphs[percentile].Draw("LP SAME");
    legend.AddEntry(&graphs[percentile], label.c_str(), "lp");
  }
  legend.Draw();
  canvas.Write();
  for (auto& graph : graphs) graph.Write();
}

}  // namespace

void accumulate_event(Outputs& outputs, const EventChannels& event_channels) {
  // Select channels within one meter and scan every MIP threshold.
  for (const auto& event_channel : event_channels) {
    if (std::hypot(event_channel.x_mm, event_channel.y_mm) > kMaxRadiusMM) continue;
    const double mip_gev = rates::mip_energy_gev(event_channel.channel.rlayerz);

    // Accumulate inclusive and origin-separated passes at each threshold.
    for (std::size_t threshold = 0; threshold < kMIPCoefficients.size(); ++threshold) {
      if (event_channel.energy_gev <= kMIPCoefficients[threshold] * mip_gev) continue;
      ++outputs.all_passes[threshold][event_channel.channel];

      // Attribute the passing channel to each tracked origin with deposited energy.
      for (std::size_t origin = 0; origin < event_channel.energy_by_origin.size(); ++origin) {
        if (event_channel.energy_by_origin[origin] > 0.0) ++outputs.origin_passes[origin][threshold][event_channel.channel];
      }
    }
  }
}

void write_output(TDirectory* parent, const Outputs& outputs, std::uint64_t n_events) {
  // Create continuous threshold scans and selected-threshold distributions.
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  auto* scan_directory = parent->mkdir("channel_hit_rate_vs_mip_r1000");
  auto* distribution_directory = parent->mkdir("channel_hit_rate_r1000");

  // Loop over the inclusive, DIS, and proton-beam-gas source groups.
  for (int source = 0; source < static_cast<int>(kSourceNames.size()); ++source) {
    auto* source_scan_directory = scan_directory->mkdir(kSourceNames[source]);
    auto* source_distribution_directory = distribution_directory->mkdir(kSourceNames[source]);
    // Write p95 and p99 threshold scans for every readout layer.
    for (int layer = 0; layer < kNReadoutLayers; ++layer) {
      write_graphs(source_scan_directory,
                   kSourceNames[source],
                   layer,
                   percentile_values(outputs, source, layer, total_time_sec));
    }

    // Write channel-rate distributions at the three comparison thresholds.
    for (std::size_t threshold : kDistributionThresholds) {
      auto* threshold_directory = source_distribution_directory->mkdir(
          ("thr" + std::to_string(static_cast<int>(100.0 * kMIPCoefficients[threshold] + 0.5)) + "MIP").c_str());

      // Build one distribution for every readout layer.
      for (int layer = 0; layer < kNReadoutLayers; ++layer) {
        TH1D histogram(("h_layer" + std::to_string(layer + 1)).c_str(),
                       ("LFHCAL " + std::string(kSourceNames[source]) + " channel rates, R < 1000 mm;rate [Hz];channels").c_str(),
                       240, 0.0, 2.0e5);
        // Fill the distribution from channels in this layer.
        for (const auto& [channel, passes] : source_map(outputs, source, threshold)) {
          if (channel.rlayerz == layer) histogram.Fill(static_cast<double>(passes) / total_time_sec);
        }
        threshold_directory->cd();
        histogram.Write();
      }
    }
  }
}

}  // namespace rates::lfhcal_analysis::occupancy::central
