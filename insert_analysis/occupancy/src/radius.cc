#include "radius.h"

#include "shared.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

namespace rates::insert_occupancy::radius {
namespace {
// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
constexpr std::array<double, 3> kCoefficients = {0.1, 0.5, 1.5};
constexpr std::array<double, 2> kPercentiles = {0.95, 0.99};
constexpr int kRadiusBins = 8;
constexpr double kRadiusMaxMm = 750.0;
const std::array<int, 2> kColors = {kBlue + 1, kRed + 1};
const std::array<int, 2> kMarkers = {20, 21};
constexpr double kMipGeV = 4e-4;

struct RadiusBinStats {
  std::vector<double> rates_gbps;
};

// ----------------------------------------------------------------------------------
// Plot helper(s)
// ----------------------------------------------------------------------------------
// Turn a percentile request into the corresponding sorted-data index.
std::size_t percentile_index(std::size_t n_values, double percentile) {
  if (n_values == 0) return 0;
  return static_cast<std::size_t>(percentile * static_cast<double>(n_values - 1));
}

std::string threshold_tag(double threshold_mip) {
  return "thr" + std::to_string(static_cast<int>(threshold_mip * 100.0 + 0.5)) + "MIP";
}

std::string percentile_tag(double percentile) {
  return "p" + std::to_string(static_cast<int>(percentile * 100.0 + 0.5));
}

std::size_t radius_bin(double radius_mm) {
  // Put everything beyond the plotting window into the last visible bin.
  if (radius_mm >= kRadiusMaxMm) return kRadiusBins - 1;
  const double fraction = radius_mm / kRadiusMaxMm;
  return std::min<std::size_t>(kRadiusBins - 1, static_cast<std::size_t>(fraction * kRadiusBins));
}

// Called when writing one radius-scan overlay canvas.
void draw_overlay(TDirectory* dir,
                  double threshold_mip,
                  const std::array<std::array<double, kRadiusBins>, 2>& values,
                  const std::array<double, kRadiusBins>& radius_centers_mm) {
  dir->cd();

  TCanvas canvas(("c_" + threshold_tag(threshold_mip) + "_data_rate_vs_radius").c_str(),
                 ("Virtual chip tail data rates vs radius, threshold " + std::to_string(threshold_mip) + " MIP;radius [mm];data rate [Gb/s]").c_str(),
                 1000,
                 800);
  canvas.SetGrid();

  auto* frame = canvas.DrawFrame(0.0, 0.0, kRadiusMaxMm, 0.6);
  frame->SetTitle(("Virtual chip tail data rates vs radius, threshold " + std::to_string(threshold_mip) + " MIP;radius [mm];data rate [Gb/s]").c_str());
  frame->SetStats(false);

  TLegend legend(0.65, 0.76, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.04);

  std::array<TGraph, 2> graphs = {
      TGraph(kRadiusBins, radius_centers_mm.data(), values[0].data()),
      TGraph(kRadiusBins, radius_centers_mm.data(), values[1].data())};

  for (int i = static_cast<int>(kPercentiles.size()) - 1; i >= 0; --i) {
    graphs[i].SetName(("g_" + threshold_tag(threshold_mip) + "_" + percentile_tag(kPercentiles[i]) + "_data_rate_vs_radius").c_str());
    graphs[i].SetLineColor(kColors[i]);
    graphs[i].SetMarkerColor(kColors[i]);
    graphs[i].SetMarkerStyle(kMarkers[i]);
    graphs[i].SetLineWidth(2);
    graphs[i].Draw("LP SAME");
    legend.AddEntry(&graphs[i], percentile_tag(kPercentiles[i]).c_str(), "lp");
  }

  legend.Draw();
  canvas.Write();
  for (auto& graph : graphs) graph.Write();
}

}  // namespace

// Called per event after the hit loop to accumulate radius-scan products.
void accumulate_event(std::array<ThresholdAccum, 3>& products,
                      const mip::EventEnergyMap& event_energy,
                      const rates::InsertToLFHCALMapper& mapper) {
  // Each threshold gets its own per-event map from virtual chip to fired-channel count.
  std::unordered_map<rates::VirtualLFHCALChipID, int, rates::VirtualLFHCALChipIDHash> active_counts[3];

  // Loop over all virtual channels in this event.
  for (const auto& [channel, energy] : event_energy) {

    // For each virtual channel, test all requested MIP thresholds.
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

      // Skip if this virtual channel does not pass the threshold.
      if (energy <= kCoefficients[threshold_index] * kMipGeV) continue;
      ++active_counts[threshold_index][mapper.chip(channel)];  // Count one passing channel on this chip for this event.
    }
  }

  // After the event is summed, convert fired virtual channels per chip into payload bits.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    for (const auto& [chip, active_count] : active_counts[threshold_index]) {
      const double event_bits =
          (kOverheadBits + kBitsPerHit * static_cast<double>(active_count)) * kSamplesPerEvent;
      products[threshold_index].bits[chip] += event_bits;
    }
  }
}

// Called once after the event loop to write radius-scan products.
void write_output(TFile& output,
                  const std::array<ThresholdAccum, 3>& products,
                  std::uint64_t n_events) {
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<double, kRadiusBins> radius_centers_mm{};
  const double width = kRadiusMaxMm / static_cast<double>(kRadiusBins);

  // Use uniform radial bins for all threshold overlays.
  for (int i = 0; i < kRadiusBins; ++i) {
    radius_centers_mm[i] = (static_cast<double>(i) + 0.5) * width;
  }

  auto* parent = output.mkdir("data_rate_vs_radius");

  // Loop over all thresholds and convert per-chip data rates into radius-binned percentiles.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    std::array<RadiusBinStats, kRadiusBins> bins{};

    // For each threshold, loop over all virtual chips.
    for (const auto& [chip, bits] : products[threshold_index].bits) {

      // Compute the chip radius and append its data rate into the corresponding radius bin.
      const double x_mm = (static_cast<double>(chip.ix) + 0.5) * 100.0;
      const double y_mm = (static_cast<double>(chip.iy) + 0.5) * 100.0;
      const double rate_gbps = (bits / total_time_sec) / 1.0e9;
      bins[radius_bin(std::hypot(x_mm, y_mm))].rates_gbps.push_back(rate_gbps);
    }

    std::array<std::array<double, kRadiusBins>, 2> values{};

    // For this threshold, reduce each radius bin into p95 and p99.
    for (int bin = 0; bin < kRadiusBins; ++bin) {
      auto& rates_gbps = bins[bin].rates_gbps;
      if (rates_gbps.empty()) continue;

      std::sort(rates_gbps.begin(), rates_gbps.end());
      values[0][bin] = rates_gbps[percentile_index(rates_gbps.size(), 0.95)];
      values[1][bin] = rates_gbps[percentile_index(rates_gbps.size(), 0.99)];
    }

    auto* dir = parent->mkdir(threshold_tag(kCoefficients[threshold_index]).c_str());
    draw_overlay(dir, kCoefficients[threshold_index], values, radius_centers_mm);
  }
}

}  // namespace rates::insert_occupancy::radius
