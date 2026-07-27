#include "occupancy/include/radius.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

namespace rates::lfhcal_analysis::occupancy::radius {
namespace {

constexpr int kRadiusBins = 30;
constexpr double kRadiusMaxMM = 3000.0;
constexpr std::array<double, 2> kPercentiles = {0.95, 0.99};

std::size_t radius_bin(double radius_mm) {
  if (radius_mm >= kRadiusMaxMM) return kRadiusBins - 1;
  return std::min<std::size_t>(kRadiusBins - 1, static_cast<std::size_t>(radius_mm / kRadiusMaxMM * kRadiusBins));
}

std::string threshold_tag(double threshold) {
  return "thr" + std::to_string(static_cast<int>(100.0 * threshold + 0.5)) + "MIP";
}

}  // namespace

void accumulate_event(std::array<ThresholdSum, kThresholds.size()>& threshold_sums,
                      const EventChannels& event_channels) {
  // Count passing channels per chip for each radius-study threshold.
  std::array<std::unordered_map<rates::LFHCALChipID, int, rates::LFHCALChipIDHash>, kThresholds.size()> active_counts;

  // Scan channels against all radius-plot thresholds.
  for (const auto& event_channel : event_channels) {
    const double mip_gev = rates::mip_energy_gev(event_channel.channel.rlayerz);
    for (std::size_t index = 0; index < kThresholds.size(); ++index) {
      if (event_channel.energy_gev <= kThresholds[index] * mip_gev) continue;
      ++active_counts[index][event_channel.chip];
    }
  }

  // Convert each fired chip into framed payload bits.
  for (std::size_t index = 0; index < kThresholds.size(); ++index) {
    for (const auto& [chip_id, active_count] : active_counts[index]) {
      threshold_sums[index].payload_bits[chip_id] +=
          (kOverheadBits + kBitsPerHit * static_cast<double>(active_count)) * kSamplesPerEvent;
    }
  }
}

void write_output(TDirectory* parent,
                  const std::array<ThresholdSum, kThresholds.size()>& threshold_sums,
                  const rates::LFHCALCellIDDecoder& decoder,
                  std::uint64_t n_events) {
  // Define common radius bins and live-time normalization.
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  auto* radius_dir = parent->mkdir("data_rate_vs_radius");
  std::array<double, kRadiusBins> centers{};
  for (int bin = 0; bin < kRadiusBins; ++bin) centers[bin] = (static_cast<double>(bin) + 0.5) * kRadiusMaxMM / kRadiusBins;

  // Build one radial rate summary for each threshold.
  for (std::size_t threshold_index = 0; threshold_index < kThresholds.size(); ++threshold_index) {
    std::array<std::vector<double>, kRadiusBins> bins;

    // Place every chip rate into its radial bin.
    for (const auto& [chip_id, payload_bits] : threshold_sums[threshold_index].payload_bits) {
      const auto position = decoder.position(chip_id);
      bins[radius_bin(std::hypot(position.x_mm, position.y_mm))].push_back(payload_bits / (total_time_sec * 1.0e9));
    }

    std::array<std::array<double, kRadiusBins>, 3> values{};

    // Compute the mean, p95, and p99 rate in every populated bin.
    for (int bin = 0; bin < kRadiusBins; ++bin) {
      auto& rates = bins[bin];
      if (rates.empty()) continue;
      double sum = 0.0;
      for (double rate : rates) sum += rate;
      values[0][bin] = sum / static_cast<double>(rates.size());
      std::sort(rates.begin(), rates.end());
      for (std::size_t percentile = 0; percentile < kPercentiles.size(); ++percentile) {
        values[percentile + 1][bin] = rates[percentile_index(rates.size(), kPercentiles[percentile])];
      }
    }

    auto* directory = radius_dir->mkdir(threshold_tag(kThresholds[threshold_index]).c_str());
    auto* graph_directory = directory->mkdir("graphs");
    directory->cd();
    const std::string title = "LFHCAL chip data rate vs radius, " + std::to_string(kThresholds[threshold_index]) +
                              " MIP;radius [mm];data rate [Gb/s]";
    TCanvas canvas(("c_" + threshold_tag(kThresholds[threshold_index]) + "_data_rate_vs_radius").c_str(), title.c_str(), 1000, 800);
    canvas.SetGrid();
    TLegend legend(0.7, 0.76, 0.88, 0.88);
    std::array<TGraph, 3> graphs = {
        TGraph(kRadiusBins, centers.data(), values[0].data()),
        TGraph(kRadiusBins, centers.data(), values[1].data()),
        TGraph(kRadiusBins, centers.data(), values[2].data())};
    TMultiGraph overlay;
        
    // Style, label, and write the three radial summary graphs.
    for (std::size_t percentile = 0; percentile < graphs.size(); ++percentile) {
      const std::string label = percentile == 0 ? "mean" : "p" + std::to_string(static_cast<int>(100.0 * kPercentiles[percentile - 1]));
      graphs[percentile].SetName(("g_" + threshold_tag(kThresholds[threshold_index]) + "_" + label).c_str());
      const int colors[3] = {kBlack, kBlue + 1, kRed + 1};
      graphs[percentile].SetLineColor(colors[percentile]);
      graphs[percentile].SetMarkerColor(graphs[percentile].GetLineColor());
      graphs[percentile].SetMarkerStyle(20 + static_cast<int>(percentile));
      overlay.Add(&graphs[percentile], "LP");
      legend.AddEntry(&graphs[percentile], label.c_str(), "lp");
    }
    overlay.SetTitle(title.c_str());
    overlay.Draw("A");
    rates::pad_axes(overlay);
    legend.Draw();
    directory->cd();
    canvas.Write();
    graph_directory->cd();
    for (auto& graph : graphs) graph.Write();
  }
}

}  // namespace rates::lfhcal_analysis::occupancy::radius
