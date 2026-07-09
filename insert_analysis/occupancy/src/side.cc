#include "side.h"

#include "shared.h"

#include <TCanvas.h>
#include <TGraph.h>

#include <array>
#include <string>
#include <unordered_map>

namespace rates::insert_occupancy::side {
namespace {
// ----------------------------------------------------------------------------------
// Constants
// ----------------------------------------------------------------------------------
constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

// Called when writing one side-level data-rate scan.
void draw_side_canvas(TDirectory* parent,
                      int side,
                      const std::array<double, kCoefficients.size()>& rates_gbps) {
  const std::string side_name = side == 0 ? "left" : "right";
  const double y_min = 0.0;
  const double y_max = 12.0;
  auto* dir = parent->mkdir(side_name.c_str());
  dir->cd();

  TCanvas canvas(("c_" + side_name + "_data_rate_vs_mip").c_str(),
                 ("Insert " + side_name + " side virtual-chip tail data rates vs MIP coefficient;MIP coefficient;data rate [Gb/s]").c_str(),
                 1000,
                 800);
  canvas.SetGrid();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), y_min, kCoefficients.back(), y_max);
  frame->SetTitle(("Insert " + side_name + " side virtual-chip tail data rates vs MIP coefficient;MIP coefficient;data rate [Gb/s]").c_str());
  frame->SetStats(false);

  TGraph graph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), rates_gbps.data());
  graph.SetName(("g_" + side_name + "_data_rate_vs_mip").c_str());
  graph.SetLineWidth(2);
  graph.SetMarkerStyle(20);
  graph.Draw("LP SAME");

  canvas.Write();
  graph.Write();
}

}  // namespace

// Called per event after the hit loop to accumulate side-level data-rate products.
void accumulate_event(std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                      const std::array<mip::EventEnergyMap, 2>& side_channel_energy_sum,
                      const rates::InsertToLFHCALMapper& mapper) {
  for (int side = 0; side < 2; ++side) {
    std::unordered_map<rates::VirtualLFHCALChipID, int, rates::VirtualLFHCALChipIDHash> active_counts[16];

    // Loop over all virtual channels on this side in this event.
    for (const auto& [channel, energy] : side_channel_energy_sum[side]) {

      // For each virtual channel, test all MIP thresholds.
      for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

        // Skip if this virtual channel does not pass the threshold.
        const int segment = rates::insert_occupancy::segment_index(channel.layer);
        if (energy <= kCoefficients[threshold_index] * rates::insert_occupancy::channel_mip_energy_gev(segment)) continue;
        ++active_counts[threshold_index][mapper.chip(channel)];
      }
    }

    // Sum the full side's virtual-chip payload for this event.
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
      double event_bits = 0.0;
      for (const auto& [chip, active_channel_count] : active_counts[threshold_index]) {
        (void)chip;
        event_bits += (kOverheadBits + kBitsPerHit * static_cast<double>(active_channel_count)) * kSamplesPerEvent;
      }
      threshold_sums[side][threshold_index].total_payload_bits += event_bits;
    }
  }
}

// Called once after the event loop to write side-level data-rate products.
void write_output(TFile& output,
                  const std::array<std::array<ThresholdSum, 16>, 2>& threshold_sums,
                  std::uint64_t n_events) {
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<std::array<double, kCoefficients.size()>, 2> rates_gbps{};

  // Convert accumulated side payloads into side data rates.
  for (int side = 0; side < 2; ++side) {
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
      rates_gbps[side][threshold_index] = (threshold_sums[side][threshold_index].total_payload_bits / total_time_sec) / 1.0e9;
    }
  }

  auto* parent = output.mkdir("side_data_rate_vs_mip");
  draw_side_canvas(parent, 0, rates_gbps[0]);
  draw_side_canvas(parent, 1, rates_gbps[1]);
}

}  // namespace rates::insert_occupancy::side
