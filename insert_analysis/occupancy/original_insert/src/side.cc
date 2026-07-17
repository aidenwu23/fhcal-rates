#include "original_insert/include/side.h"

#include "original_insert/include/shared.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

#include <array>
#include <unordered_map>

namespace rates::insert_analysis::original_insert::side {
namespace {

constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
constexpr double kYMaxGbps = 25.0;

void draw_side_canvas(TDirectory* dir,
                      const std::array<std::array<double, kCoefficients.size()>, 2>& rates_gbps) {
  dir->cd();
  TCanvas canvas("c_side_data_rate_vs_mip",
                 "Insert total data rate by side vs MIP coefficient (original insert);MIP coefficient;data rate [Gb/s]",
                 1000,
                 800);
  canvas.SetGrid();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), 0.0, kCoefficients.back(), kYMaxGbps);
  frame->SetTitle("Insert total data rate by side vs MIP coefficient (original insert);MIP coefficient;data rate [Gb/s]");
  frame->GetXaxis()->SetTitle("MIP coefficient");
  frame->GetYaxis()->SetTitle("data rate [Gb/s]");
  frame->SetStats(false);

  TLegend legend(0.65, 0.76, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.04);

  TGraph left_graph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), rates_gbps[0].data());
  left_graph.SetName("g_left_data_rate_vs_mip");
  left_graph.SetLineWidth(2);
  left_graph.SetLineColor(kBlue + 1);
  left_graph.SetMarkerColor(kBlue + 1);
  left_graph.SetMarkerStyle(20);
  left_graph.Draw("LP SAME");
  legend.AddEntry(&left_graph, "left", "lp");

  TGraph right_graph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), rates_gbps[1].data());
  right_graph.SetName("g_right_data_rate_vs_mip");
  right_graph.SetLineWidth(2);
  right_graph.SetLineColor(kRed + 1);
  right_graph.SetMarkerColor(kRed + 1);
  right_graph.SetMarkerStyle(21);
  right_graph.Draw("LP SAME");
  legend.AddEntry(&right_graph, "right", "lp");

  legend.Draw();
  canvas.Write();
  left_graph.Write();
  right_graph.Write();
}

}  // namespace

void accumulate_event(std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                      const OccupancyMode& mode,
                      const std::array<mip::EventEnergyMap, 2>& side_channel_energies,
                      const rates::InsertToLFHCALMapper& mapper) {

  // Do both sides.
  for (int side = 0; side < 2; ++side) {
    std::unordered_map<rates::VirtualLFHCALPizzaChipID,
                       int,
                       rates::VirtualLFHCALPizzaChipIDHash>
        active_counts[16];

    // Per side, loop through all channels.
    for (const auto& [channel, energy] : side_channel_energies[side]) {

      // Per channel, test every threshold.
      for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
        if (energy <= kCoefficients[threshold_index] * channel_mip_energy_gev(mode, channel.layer)) continue;
        ++active_counts[threshold_index][mapper.pizza_chip(channel.layer, channel.x_mm, channel.y_mm, side)];
      }
    }

    // Per side, loop over all thresholds.
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
      double payload_bits = 0.0;

      // Per threshold, sum the payload from every active chip.
      for (const auto& [chip, active_channel_count] : active_counts[threshold_index]) {
        (void)chip;
        payload_bits +=
            (kOverheadBits + kBitsPerHit * static_cast<double>(active_channel_count)) * kSamplesPerEvent;
      }
      threshold_sums[side][threshold_index].total_payload_bits += payload_bits;
    }
  }
}

void write_output(TDirectory* parent,
                  const std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                  std::uint64_t n_events) {
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<std::array<double, kCoefficients.size()>, 2> rates_gbps{};

  // Loop over both sides.
  for (int side = 0; side < 2; ++side) {

    // Per side loop over all thresholds.
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

      // Compute rate.
      rates_gbps[side][threshold_index] = (threshold_sums[side][threshold_index].total_payload_bits / total_time_sec) / 1.0e9;
    }
  }

  auto* dir = parent->mkdir("side_data_rate_vs_mip");
  draw_side_canvas(dir, rates_gbps);
}

}  // namespace rates::insert_analysis::original_insert::side
