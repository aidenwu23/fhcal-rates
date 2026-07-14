#include "full_lfhcal/include/mip.h"

#include "full_lfhcal/include/shared.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

#include <algorithm>
#include <array>
#include <string>
#include <unordered_set>
#include <vector>

namespace rates::insert_analysis::full_lfhcal::mip {
namespace {
// ----------------------------------------------------------------------------------
// Constants
// ----------------------------------------------------------------------------------
constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
constexpr std::array<double, 2> kPercentiles = {0.95, 0.99};
const std::array<int, 2> kColors = {kBlue + 1, kRed + 1};
const std::array<int, 2> kMarkers = {20, 21};

// ----------------------------------------------------------------------------------
// Plot helper(s)
// ----------------------------------------------------------------------------------
std::string percentile_tag(double percentile) {
  return "p" + std::to_string(static_cast<int>(percentile * 100.0 + 0.5));
}

// Called when writing one MIP-scan overlay canvas.
void draw_overlay(TDirectory* dir,
                  const char* canvas_name,
                  const char* title,
                  const char* graph_prefix,
                  const char* y_title,
                  double y_min,
                  double y_max,
                  const std::array<std::array<double, 16>, 2>& percentile_values) {
  dir->cd();

  TCanvas canvas(canvas_name, title, 1000, 800);
  canvas.SetGrid();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), y_min, kCoefficients.back(), y_max);
  frame->SetTitle(title);
  frame->GetXaxis()->SetTitle("MIP coefficient");
  frame->GetYaxis()->SetTitle(y_title);
  frame->SetStats(false);

  TLegend legend(0.65, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.04);

  std::array<TGraph, 2> graphs = {
      TGraph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), percentile_values[0].data()),
      TGraph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), percentile_values[1].data())};

  // Draw one curve per requested percentile.
  for (int i = static_cast<int>(kPercentiles.size()) - 1; i >= 0; --i) {
    graphs[i].SetName((std::string(graph_prefix) + "_" + percentile_tag(kPercentiles[i])).c_str());
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

// Called per event after the hit loop to accumulate MIP-threshold products.
void accumulate_event(std::array<ChannelThresholdSum, 16>& channel_threshold_sums,
                      std::array<ChannelDataThresholdSum, 16>& channel_data_threshold_sums,
                      std::array<ChipThresholdSum, 16>& chip_threshold_sums,
                      std::array<ChipDataThresholdSum, 16>& chip_data_threshold_sums,
                      const OccupancyMode& mode,
                      const EventEnergyMap& channel_energies,
                      const rates::InsertToLFHCALMapper& mapper) {
  // These per-threshold containers exist for one event, then fold into the full sample.
  std::unordered_map<rates::VirtualLFHCALChipID,
                     std::unordered_set<rates::VirtualLFHCALChannelID, rates::VirtualLFHCALChannelIDHash>,
                     rates::VirtualLFHCALChipIDHash>
      fired_channels_by_chip[16];
  std::unordered_set<rates::VirtualLFHCALChipID, rates::VirtualLFHCALChipIDHash> fired_chips[16];

  // Loop over all virtual channels in this event.
  for (const auto& [channel, energy] : channel_energies) {

    // For each virtual channel, test all MIP thresholds.
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

      // Skip if this virtual channel does not pass the threshold.
      if (energy <= kCoefficients[threshold_index] * rates::insert_analysis::full_lfhcal::channel_mip_energy_gev(mode, channel.layer)) continue;

      // Else, count one passing event and its channel payload bits.
      ++channel_threshold_sums[threshold_index].pass_counts[channel];
      channel_data_threshold_sums[threshold_index].payload_bits[channel] += kBitsPerHit * kSamplesPerEvent;

      // Map the channel to a chip.
      const auto chip = mapper.chip(channel);
      fired_chips[threshold_index].insert(chip); // Track that some channel on this chip fired at least once per event.
      fired_channels_by_chip[threshold_index][chip].insert(channel);  // Track how many distinct channels fired on that chip.
    }
  }

  // After the event is summed, count fired chips and accumulate chip payload bits.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    for (const auto& chip : fired_chips[threshold_index]) {
      ++chip_threshold_sums[threshold_index].pass_counts[chip];
    }
    for (const auto& [chip, fired_channels] : fired_channels_by_chip[threshold_index]) {
      const double event_bits =
          (kOverheadBits + kBitsPerHit * static_cast<double>(fired_channels.size())) * kSamplesPerEvent;
      chip_data_threshold_sums[threshold_index].payload_bits[chip] += event_bits;
    }
  }
}

// Called once after the event loop to write MIP-scan products.
void write_output(TDirectory* parent,
                  const std::array<ChannelThresholdSum, 16>& channel_threshold_sums,
                  const std::array<ChannelDataThresholdSum, 16>& channel_data_threshold_sums,
                  const std::array<ChipThresholdSum, 16>& chip_threshold_sums,
                  const std::array<ChipDataThresholdSum, 16>& chip_data_threshold_sums,
                  const OccupancyMode& mode,
                  std::uint64_t n_events) {
  // Convert full-sample counts into rates using the simulated event window.
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<std::array<double, 16>, 2> channel_hit_rate_percentile_values{};
  std::array<std::array<double, 16>, 2> channel_data_rate_percentile_values{};
  std::array<std::array<double, 16>, 2> chip_hit_rate_percentile_values{};
  std::array<std::array<double, 16>, 2> chip_data_rate_percentile_values{};

  // Loop over all thresholds and convert accumulated counts into percentile curves.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    std::vector<double> channel_rates_hz;

    // Collect one virtual-channel rate for every channel that ever passed.
    for (const auto& [channel, pass_count] : channel_threshold_sums[threshold_index].pass_counts) {
      (void)channel;
      channel_rates_hz.push_back(static_cast<double>(pass_count) / total_time_sec);
    }
    if (!channel_rates_hz.empty()) {
      std::sort(channel_rates_hz.begin(), channel_rates_hz.end());
      for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
        channel_hit_rate_percentile_values[percentile_slot][threshold_index] =
            channel_rates_hz[rates::insert_analysis::full_lfhcal::percentile_index(channel_rates_hz.size(), kPercentiles[percentile_slot])];
      }
    }

    // Build the matching data-rate distribution from accumulated channel payload bits.
    std::vector<double> channel_data_rates_gbps;
    for (const auto& [channel, payload_bits] : channel_data_threshold_sums[threshold_index].payload_bits) {
      (void)channel;
      channel_data_rates_gbps.push_back((payload_bits / total_time_sec) / 1.0e9);
    }
    if (!channel_data_rates_gbps.empty()) {
      std::sort(channel_data_rates_gbps.begin(), channel_data_rates_gbps.end());
      for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
        channel_data_rate_percentile_values[percentile_slot][threshold_index] =
            channel_data_rates_gbps[rates::insert_analysis::full_lfhcal::percentile_index(channel_data_rates_gbps.size(), kPercentiles[percentile_slot])];
      }
    }

    // Repeat the percentile reduction for virtual chips.
    std::vector<double> chip_rates_hz;

    // Collect one virtual-chip hit rate for every chip that ever passed.
    for (const auto& [chip, pass_count] : chip_threshold_sums[threshold_index].pass_counts) {
      (void)chip;
      chip_rates_hz.push_back(static_cast<double>(pass_count) / total_time_sec);
    }
    if (!chip_rates_hz.empty()) {
      std::sort(chip_rates_hz.begin(), chip_rates_hz.end());

      // Do it for all requested percentiles.
      for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
        chip_hit_rate_percentile_values[percentile_slot][threshold_index] =
            chip_rates_hz[rates::insert_analysis::full_lfhcal::percentile_index(chip_rates_hz.size(), kPercentiles[percentile_slot])];
      }
    }

    std::vector<double> chip_data_rates_gbps;

    for (const auto& [chip, payload_bits] : chip_data_threshold_sums[threshold_index].payload_bits) {
      (void)chip;
      chip_data_rates_gbps.push_back((payload_bits / total_time_sec) / 1.0e9);
    }
    if (!chip_data_rates_gbps.empty()) {
      std::sort(chip_data_rates_gbps.begin(), chip_data_rates_gbps.end());
      for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
        chip_data_rate_percentile_values[percentile_slot][threshold_index] =
            chip_data_rates_gbps[rates::insert_analysis::full_lfhcal::percentile_index(chip_data_rates_gbps.size(), kPercentiles[percentile_slot])];
      }
    }
  }

  // Write separate percentile overlays for channel and chip hit and data rates.
  auto* channel_hit_dir = parent->mkdir("channel_hit_rate_vs_mip");
  draw_overlay(channel_hit_dir,
               "c_channel_hit_rate_vs_mip",
               "Insert channel tail hit rates vs MIP coefficient (full LFHCal readout);MIP coefficient;rate [Hz]",
               "g_channel_hit_rate_vs_mip",
               "rate [Hz]",
               0.0,
               2.0e5,
               channel_hit_rate_percentile_values);

  auto* channel_data_dir = parent->mkdir("channel_data_rate_vs_mip");
  draw_overlay(channel_data_dir,
               "c_channel_data_rate_vs_mip",
               "Insert channel tail data rates vs MIP coefficient (full LFHCal readout);MIP coefficient;data rate [Gb/s]",
               "g_channel_data_rate_vs_mip",
               "data rate [Gb/s]",
               0.0,
               mode.drop_inner_ring ? 25.5e-3 : 26.0e-3,
               channel_data_rate_percentile_values);

  auto* chip_hit_dir = parent->mkdir("chip_hit_rate_vs_mip");
  draw_overlay(chip_hit_dir,
               "c_chip_hit_rate_vs_mip",
               "Insert chip tail hit rates vs MIP coefficient (full LFHCal readout);MIP coefficient;rate [Hz]",
               "g_chip_hit_rate_vs_mip",
               "rate [Hz]",
               8.0e4,
               4.2e5,
               chip_hit_rate_percentile_values);

  auto* chip_data_dir = parent->mkdir("chip_data_rate_vs_mip");
  draw_overlay(chip_data_dir,
               "c_chip_data_rate_vs_mip",
               "Insert chip tail data rates vs MIP coefficient (full LFHCal readout);MIP coefficient;data rate [Gb/s]",
               "g_chip_data_rate_vs_mip",
               "data rate [Gb/s]",
               0.0,
               0.75,
               chip_data_rate_percentile_values);
}

}  // namespace rates::insert_analysis::full_lfhcal::mip
