#include "lfhcal_segments/include/side.h"

#include "lfhcal_segments/include/shared.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

#include <array>
#include <string>
#include <unordered_map>

namespace rates::insert_analysis::lfhcal_segments::side {
namespace {
// ----------------------------------------------------------------------------------
// Constants
// ----------------------------------------------------------------------------------
constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

// Called when writing the side-level data-rate scan.
void draw_side_canvas(TDirectory* dir,
                      const std::array<std::array<double, kCoefficients.size()>, 2>& rates_gbps,
                      double y_max) {
  const double y_min = 0.0;
  dir->cd();

  TCanvas canvas("c_side_data_rate_vs_mip",
                 "Insert side virtual-chip tail data rates vs MIP coefficient;MIP coefficient;data rate [Gb/s]",
                 1000,
                 800);
  canvas.SetGrid();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), y_min, kCoefficients.back(), y_max);
  frame->SetTitle("Insert side virtual-chip tail data rates vs MIP coefficient;MIP coefficient;data rate [Gb/s]");
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

// Called per event after the hit loop to accumulate side-level data-rate products.
void accumulate_event(std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                      const OccupancyMode& mode,
                      const std::array<mip::EventEnergyMap, 2>& side_channel_energies,
                      const rates::InsertToLFHCALMapper& mapper) {
  (void)mode;
  for (int side = 0; side < 2; ++side) {
    std::unordered_map<rates::VirtualLFHCALChipID, int, rates::VirtualLFHCALChipIDHash> active_counts[16];

    // Loop over all virtual channels on this side in this event.
    for (const auto& [channel, energy] : side_channel_energies[side]) {

      // For each virtual channel, test all MIP thresholds.
      for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

        // Skip if this virtual channel does not pass the threshold.
        if (energy <= kCoefficients[threshold_index] * rates::insert_analysis::lfhcal_segments::channel_mip_energy_gev(mode, channel.layer)) continue;
        ++active_counts[threshold_index][mapper.chip(channel)];
      }
    }

    // Sum the full side's virtual-chip payload for this event.
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
      double payload_bits = 0.0;
      for (const auto& [chip, active_channel_count] : active_counts[threshold_index]) {
        (void)chip;
        payload_bits += (kOverheadBits + kBitsPerHit * static_cast<double>(active_channel_count)) * kSamplesPerEvent;
      }
      threshold_sums[side][threshold_index].total_payload_bits += payload_bits; // Add this event's side payload.
    }
  }
}

// Called once after the event loop to write side-level data-rate products.
void write_output(TDirectory* parent,
                  const std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                  const OccupancyMode& mode,
                  std::uint64_t n_events) {
  // Convert the summed side payloads into rates using the simulated event window.
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<std::array<double, kCoefficients.size()>, 2> rates_gbps{};

  // Convert accumulated side payloads into side data rates.
  for (int side = 0; side < 2; ++side) {
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
      rates_gbps[side][threshold_index] = (threshold_sums[side][threshold_index].total_payload_bits / total_time_sec) / 1.0e9;
    }
  }

  auto* dir = parent->mkdir("side_data_rate_vs_mip");
  const double y_max = mode.drop_inner_ring ? 11.0 : 12.0;
  draw_side_canvas(dir, rates_gbps, y_max);
}

}  // namespace rates::insert_analysis::lfhcal_segments::side
