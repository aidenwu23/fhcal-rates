#include "mip.h"

#include "shared.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

#include <algorithm>
#include <array>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace rates::insert_occupancy::mip {
namespace {
// ----------------------------------------------------------------------------------
// Constants
// ----------------------------------------------------------------------------------
constexpr double kMipGeV = 4e-4;
constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
constexpr std::array<double, 2> kPercentiles = {0.95, 0.99};
const std::array<int, 2> kColors = {kBlue + 1, kRed + 1};
const std::array<int, 2> kMarkers = {20, 21};

// ----------------------------------------------------------------------------------
// Plot helper(s)
// ----------------------------------------------------------------------------------
// Turn a percentile request into the corresponding sorted-data index.
std::size_t percentile_index(std::size_t n_values, double percentile) {
  if (n_values == 0) return 0;
  return static_cast<std::size_t>(percentile * static_cast<double>(n_values - 1));
}

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
                  const std::array<std::array<double, 16>, 2>& values) {
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
      TGraph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), values[0].data()),
      TGraph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), values[1].data())};

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
void accumulate_event(std::array<ChannelThresholdAccum, 16>& channel_products,
                      std::array<ChipThresholdAccum, 16>& chip_products,
                      std::array<DataThresholdAccum, 16>& data_products,
                      const EventEnergyMap& event_energy,
                      const rates::InsertToLFHCALMapper& mapper) {
  // These per-threshold containers live for one event, then get folded into the full sample.
  std::unordered_map<rates::VirtualLFHCALChipID,
                     std::unordered_set<rates::VirtualLFHCALChannelID, rates::VirtualLFHCALChannelIDHash>,
                     rates::VirtualLFHCALChipIDHash>
      fired_channels_by_chip[16];
  std::unordered_set<rates::VirtualLFHCALChipID, rates::VirtualLFHCALChipIDHash> fired_chips[16];

  // Loop over all virtual channels in this event.
  for (const auto& [channel, energy] : event_energy) {

    // For each virtual channel, test all MIP thresholds.
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

      // Skip if this virtual channel does not pass the threshold.
      if (energy <= kCoefficients[threshold_index] * kMipGeV) continue;

      // Promote the passing channel into the full-sample counters.
      ++channel_products[threshold_index].passes[channel];
      const auto chip = mapper.chip(channel);
      fired_chips[threshold_index].insert(chip);                   // Track that this chip fired at least once this event.
      fired_channels_by_chip[threshold_index][chip].insert(channel);  // Track how many distinct channels fired on that chip.
    }
  }

  // After the event is summed, count fired chips and accumulate chip payload bits.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    for (const auto& chip : fired_chips[threshold_index]) {
      ++chip_products[threshold_index].passes[chip];
    }
    for (const auto& [chip, fired_channels] : fired_channels_by_chip[threshold_index]) {
      const double event_bits =
          (kOverheadBits + kBitsPerHit * static_cast<double>(fired_channels.size())) * kSamplesPerEvent;
      data_products[threshold_index].bits[chip] += event_bits;
    }
  }
}

// Called once after the event loop to write MIP-scan products.
void write_output(TFile& output,
                  const std::array<ChannelThresholdAccum, 16>& channel_products,
                  const std::array<ChipThresholdAccum, 16>& chip_products,
                  const std::array<DataThresholdAccum, 16>& data_products,
                  std::uint64_t n_events) {
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<std::array<double, 16>, 2> channel_values{};
  std::array<std::array<double, 16>, 2> chip_values{};
  std::array<std::array<double, 16>, 2> data_values{};

  // Loop over all thresholds and convert accumulated counts into percentile curves.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    std::vector<double> channel_rates_hz;

    // Collect one virtual-channel rate for every channel that ever passed.
    for (const auto& [channel, passes] : channel_products[threshold_index].passes) {
      (void)channel;
      channel_rates_hz.push_back(static_cast<double>(passes) / total_time_sec);
    }
    if (!channel_rates_hz.empty()) {
      std::sort(channel_rates_hz.begin(), channel_rates_hz.end());
      for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
        channel_values[percentile_slot][threshold_index] =
            channel_rates_hz[percentile_index(channel_rates_hz.size(), kPercentiles[percentile_slot])];
      }
    }

    std::vector<double> chip_rates_hz;

    // Collect one virtual-chip hit rate for every chip that ever passed.
    for (const auto& [chip, passes] : chip_products[threshold_index].passes) {
      (void)chip;
      chip_rates_hz.push_back(static_cast<double>(passes) / total_time_sec);
    }
    if (!chip_rates_hz.empty()) {
      std::sort(chip_rates_hz.begin(), chip_rates_hz.end());

      // Do it for all requested percentiles.
      for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
        chip_values[percentile_slot][threshold_index] =
            chip_rates_hz[percentile_index(chip_rates_hz.size(), kPercentiles[percentile_slot])];
      }
    }

    std::vector<double> data_rates_gbps;

    // Collect one virtual-chip data rate for every chip that ever passed.
    for (const auto& [chip, bits] : data_products[threshold_index].bits) {
      (void)chip;
      data_rates_gbps.push_back((bits / total_time_sec) / 1.0e9);
    }
    if (!data_rates_gbps.empty()) {
      std::sort(data_rates_gbps.begin(), data_rates_gbps.end());

      // Do it for all percentiles requested.
      for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
        data_values[percentile_slot][threshold_index] =
            data_rates_gbps[percentile_index(data_rates_gbps.size(), kPercentiles[percentile_slot])];
      }
    }
  }

  // Write the three scan families into separate directories in the output file.
  auto* channel_dir = output.mkdir("channel_rate_vs_mip");
  draw_overlay(channel_dir,
               "c_channel_rate_vs_mip",
               "Virtual channel rate vs MIP coefficient;MIP coefficient;rate [Hz/virtual channel]",
               "g_channel_rate_vs_mip",
               "rate [Hz/virtual channel]",
               4.0e4,
               2.0e5,
               channel_values);

  auto* chip_dir = output.mkdir("chip_rate_vs_mip");
  draw_overlay(chip_dir,
               "c_chip_rate_vs_mip",
               "Virtual chip rate vs MIP coefficient;MIP coefficient;rate [Hz/virtual chip]",
               "g_chip_rate_vs_mip",
               "rate [Hz/virtual chip]",
               2.6e5,
               4.2e5,
               chip_values);

  auto* data_dir = output.mkdir("data_rate_vs_mip");
  draw_overlay(data_dir,
               "c_data_rate_vs_mip",
               "Virtual chip tail data rates vs MIP coefficient;MIP coefficient;data rate [Gb/s]",
               "g_data_rate_vs_mip",
               "data rate [Gb/s]",
               0.2,
               0.9,
               data_values);
}

}  // namespace rates::insert_occupancy::mip
