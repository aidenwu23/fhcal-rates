#include "contour.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>

#include <string>
#include <vector>

namespace rates::insert_analysis::contour {
namespace {

// One fired chip event contains fixed framing overhead plus one payload word per active channel.
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;
constexpr int kContourColor = 17;

// Solve the chip data-rate formula for the mean number of fired channels.
double contour_y(double rate_hz, double rate_gbps) {
  return ((rate_gbps * 1.0e9 / rate_hz) / kSamplesPerEvent - kOverheadBits) / kBitsPerHit;
}

}  // namespace

void write_output(TDirectory* parent,
                  const std::vector<double>& chip_rates_hz,
                  const std::vector<double>& mean_fired_channels,
                  double x_max,
                  double y_min,
                  double y_max,
                  const std::vector<double>& contour_rates_gbps,
                  const char* title_suffix) {
  parent->cd();

  // Create one frame shared by the chip points and constant-data-rate contour lines.
  const std::string title =
      std::string("Mean fired channels per active event vs rate for each chip (") + title_suffix +
      ");rate [Hz];mean fired channels per active event";
  TCanvas canvas("c_fired_channels_vs_rate", title.c_str(), 1000, 800);
  canvas.SetGrid();

  auto* frame = canvas.DrawFrame(0.0, y_min, x_max, y_max);
  frame->SetTitle(title.c_str());
  frame->SetStats(false);

  // Each point represents one chip's firing rate and its mean channel count when active.
  TGraph graph(static_cast<int>(chip_rates_hz.size()), chip_rates_hz.data(), mean_fired_channels.data());
  graph.SetName("g_fired_channels_vs_rate");
  graph.SetMarkerStyle(20);
  graph.Draw("P SAME");

  std::vector<TGraph> contour_graphs;
  contour_graphs.reserve(contour_rates_gbps.size());

  // Build one curve for each requested total chip data rate.
  for (std::size_t contour_index = 0; contour_index < contour_rates_gbps.size(); ++contour_index) {
    std::vector<double> x_values;
    std::vector<double> y_values;

    // Step across the visible chip-rate range and solve for the corresponding channel count.
    for (int step = 1; step <= 2000; ++step) {
      const double rate_hz = x_max * static_cast<double>(step) / 2000.0;
      const double y_value = contour_y(rate_hz, contour_rates_gbps[contour_index]);
      // Retain curve points that fall inside the displayed frame.
      if (y_value < y_min || y_value > y_max) continue;
      x_values.push_back(rate_hz);
      y_values.push_back(y_value);
    }
    if (x_values.empty()) continue;

    // Draw the sampled constant-rate curve over the chip scatter plot.
    contour_graphs.emplace_back(static_cast<int>(x_values.size()), x_values.data(), y_values.data());
    auto& contour = contour_graphs.back();
    contour.SetName(("g_contour_" + std::to_string(contour_index)).c_str());
    contour.SetLineColor(kContourColor);
    contour.SetLineStyle(2);
    contour.SetLineWidth(2);
    contour.Draw("L SAME");
  }

  // Store the combined canvas and its individual graphs for later inspection.
  canvas.Write();
  graph.Write();
  for (auto& contour : contour_graphs) contour.Write();
}

}  // namespace rates::insert_analysis::contour
