#include "contour.h"

#include "occupancy/include/shared.h"

#include <TCanvas.h>
#include <TGraph.h>

#include <array>
#include <string>
#include <vector>

namespace rates::lfhcal_analysis::occupancy::contour {
namespace {

constexpr std::array<double, 7> kRatesGbps = {0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035};

double contour_y(double rate_hz, double rate_gbps) {
  return ((rate_gbps * 1.0e9 / rate_hz) / kSamplesPerEvent - kOverheadBits) / kBitsPerHit;
}

}  // namespace

void write_output(TDirectory* parent,
                  const std::vector<double>& chip_rates_hz,
                  const std::vector<double>& mean_fired_channels) {
  // Draw chip rate versus fired-channel multiplicity.
  parent->cd();
  const char* title = "LFHCAL mean fired channels per active chip-event;rate [Hz];mean fired channels";
  TCanvas canvas("c_fired_channels_vs_rate", title, 1000, 800);
  canvas.SetGrid();
  auto* frame = canvas.DrawFrame(0.0, 0.0, 5.0e4, 5.0);
  frame->SetTitle(title);

  TGraph graph(static_cast<int>(chip_rates_hz.size()), chip_rates_hz.data(), mean_fired_channels.data());
  graph.SetName("g_fired_channels_vs_rate");
  graph.SetMarkerStyle(20);
  graph.Draw("P SAME");

  std::vector<TGraph> contours;
  contours.reserve(kRatesGbps.size());
  // Overlay one constant-data-rate contour for each target bandwidth.
  for (std::size_t index = 0; index < kRatesGbps.size(); ++index) {
    std::vector<double> x;
    std::vector<double> y;
    
    // Sample the contour across the visible chip-rate range.
    for (int step = 1; step <= 400; ++step) {
      const double rate_hz = 5.0e4 * static_cast<double>(step) / 400.0;
      const double channels = contour_y(rate_hz, kRatesGbps[index]);
      if (channels < 0.0 || channels > 5.0) continue;
      x.push_back(rate_hz);
      y.push_back(channels);
    }
    if (x.empty()) continue;
    contours.emplace_back(static_cast<int>(x.size()), x.data(), y.data());
    auto& contour = contours.back();
    contour.SetName(("g_contour_" + std::to_string(index)).c_str());
    contour.SetLineColor(17);
    contour.SetLineStyle(2);
    contour.Draw("L SAME");
  }

  // Write the combined canvas and its component graphs.
  canvas.Write();
  graph.Write();
  for (auto& contour : contours) contour.Write();
}

}  // namespace rates::lfhcal_analysis::occupancy::contour
