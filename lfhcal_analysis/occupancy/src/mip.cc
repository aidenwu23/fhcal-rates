#include "occupancy/include/mip.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TLegend.h>

#include <algorithm>
#include <array>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace rates::lfhcal_analysis::occupancy::mip {
namespace {

constexpr std::array<double, 2> kPercentiles = {0.95, 0.99};
constexpr std::array<int, 2> kColors = {kBlue + 1, kRed + 1};

void draw_overlay(TDirectory* directory,
                  const char* canvas_name,
                  const char* title,
                  const char* graph_name,
                  const std::array<std::array<double, kMIPCoefficients.size()>, 2>& values) {
  // Draw p95 and p99 threshold curves on one canvas.
  directory->cd();
  TCanvas canvas(canvas_name, title, 1000, 800);
  canvas.SetGrid();
  auto* frame = canvas.DrawFrame(0.0, 0.0, 1.5, 1.1 * std::max(values[1][0], 1.0));
  frame->SetTitle(title);

  TLegend legend(0.68, 0.75, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  std::array<TGraph, 2> graphs = {
      TGraph(static_cast<int>(kMIPCoefficients.size()), kMIPCoefficients.data(), values[0].data()),
      TGraph(static_cast<int>(kMIPCoefficients.size()), kMIPCoefficients.data(), values[1].data())};

  // Style, label, and write each percentile curve.
  for (std::size_t index = 0; index < graphs.size(); ++index) {
    graphs[index].SetName((std::string(graph_name) + "_p" + std::to_string(static_cast<int>(100.0 * kPercentiles[index]))).c_str());
    graphs[index].SetLineColor(kColors[index]);
    graphs[index].SetMarkerColor(kColors[index]);
    graphs[index].SetMarkerStyle(20 + static_cast<int>(index));
    graphs[index].Draw("LP SAME");
    legend.AddEntry(&graphs[index], ("p" + std::to_string(static_cast<int>(100.0 * kPercentiles[index]))).c_str(), "lp");
  }
  legend.Draw();
  canvas.Write();
  for (auto& graph : graphs) graph.Write();
}

void write_mean_graph(TDirectory* directory,
                      const std::array<double, kMIPCoefficients.size()>& values) {
  directory->cd();
  TGraph graph(static_cast<int>(kMIPCoefficients.size()), kMIPCoefficients.data(), values.data());
  graph.SetName("g_mean_chip_data_rate_vs_mip");
  graph.SetMarkerStyle(20);
  graph.Write();
}

void write_chip_channel_graphs(TDirectory* directory,
                               const std::array<double, kMIPCoefficients.size()>& rates_hz,
                               const std::array<double, kMIPCoefficients.size()>& means) {
  directory->cd();
  TGraph rate_graph(static_cast<int>(kMIPCoefficients.size()), kMIPCoefficients.data(), rates_hz.data());
  rate_graph.SetName("g_chip_channel_rate_vs_mip");
  rate_graph.SetMarkerStyle(20);
  rate_graph.Write();
  TGraph mean_graph(static_cast<int>(kMIPCoefficients.size()), kMIPCoefficients.data(), means.data());
  mean_graph.SetName("g_active_channels_vs_mip");
  mean_graph.SetMarkerStyle(20);
  mean_graph.Write();
}

void write_total_rate_graph(TDirectory* directory,
                            const std::array<double, kMIPCoefficients.size()>& rates_gbps) {
  directory->cd();
  TGraph graph(static_cast<int>(kMIPCoefficients.size()), kMIPCoefficients.data(), rates_gbps.data());
  graph.SetName("g_lfhcal_rate");
  graph.SetMarkerStyle(20);
  graph.Write();
}

template <typename Map>
std::array<double, 2> rate_percentiles(const Map& values, double scale) {
  // Convert an accumulated map into a sorted rate distribution.
  std::vector<double> rates;
  rates.reserve(values.size());
  for (const auto& [id, value] : values) {
    (void)id;
    rates.push_back(static_cast<double>(value) * scale);
  }
  if (rates.empty()) return {};

  // Return the p95 and p99 entries from the sorted distribution.
  std::sort(rates.begin(), rates.end());
  return {
      rates[percentile_index(rates.size(), kPercentiles[0])],
      rates[percentile_index(rates.size(), kPercentiles[1])]};
}

}  // namespace

void accumulate_event(std::array<ThresholdSum, kMIPCoefficients.size()>& threshold_sums,
                      const EventChannels& event_channels) {
  // Prepare per-threshold chip channel counts and the unthresholded chip set.
  std::array<std::unordered_map<rates::LFHCALChipID, int, rates::LFHCALChipIDHash>, kMIPCoefficients.size()> active_counts;
  std::unordered_set<rates::LFHCALChipID, rates::LFHCALChipIDHash> event_chips;

  // Test every channel against every MIP threshold.
  for (const auto& event_channel : event_channels) {
    event_chips.insert(event_channel.chip);
    const double mip_gev = rates::mip_energy_gev(event_channel.channel.rlayerz);

    // Add passing channels to the matching channel and chip accumulators.
    for (std::size_t index = 0; index < kMIPCoefficients.size(); ++index) {
      if (event_channel.energy_gev <= kMIPCoefficients[index] * mip_gev) continue;
      ++threshold_sums[index].channel_passes[event_channel.channel];
      ++active_counts[index][event_channel.chip];
    }
  }

  // Convert per-event active-channel counts into chip passes and payload bits.
  for (std::size_t index = 0; index < kMIPCoefficients.size(); ++index) {

    // Accumulate every chip that passes this threshold.
    for (const auto& [chip_id, active_count] : active_counts[index]) {
      ++threshold_sums[index].chip_passes[chip_id];
      threshold_sums[index].total_active_channels += static_cast<std::uint64_t>(active_count);
      threshold_sums[index].chip_payload_bits[chip_id] +=
          (kOverheadBits + kBitsPerHit * static_cast<double>(active_count)) * kSamplesPerEvent;
    }
    threshold_sums[index].total_chip_instances += event_chips.size();
  }
}

void write_output(TDirectory* parent,
                  const std::array<ThresholdSum, kMIPCoefficients.size()>& threshold_sums,
                  std::uint64_t n_events) {
  // Allocate output arrays for percentile, mean, and total threshold scans.
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<std::array<double, kMIPCoefficients.size()>, 2> channel_hit_values{};
  std::array<std::array<double, kMIPCoefficients.size()>, 2> chip_hit_values{};
  std::array<std::array<double, kMIPCoefficients.size()>, 2> chip_data_values{};
  std::array<double, kMIPCoefficients.size()> mean_chip_data_values{};
  std::array<double, kMIPCoefficients.size()> chip_channel_rates_hz{};
  std::array<double, kMIPCoefficients.size()> mean_active_channels{};
  std::array<double, kMIPCoefficients.size()> total_data_rates_gbps{};

  // Reduce each threshold's accumulated maps into plotted values.
  for (std::size_t index = 0; index < kMIPCoefficients.size(); ++index) {
    const auto channel_rates = rate_percentiles(threshold_sums[index].channel_passes, 1.0 / total_time_sec);
    const auto chip_rates = rate_percentiles(threshold_sums[index].chip_passes, 1.0 / total_time_sec);
    const auto chip_data_rates = rate_percentiles(threshold_sums[index].chip_payload_bits, 1.0 / (total_time_sec * 1.0e9));

    // Compute mean per-chip and total-detector data rates.
    if (!threshold_sums[index].chip_payload_bits.empty()) {
      double total_payload_bits = 0.0;

      // Sum payload across every chip seen at this threshold.
      for (const auto& [chip_id, payload_bits] : threshold_sums[index].chip_payload_bits) {
        (void)chip_id;
        total_payload_bits += payload_bits;
      }
      mean_chip_data_values[index] = total_payload_bits /
          (static_cast<double>(threshold_sums[index].chip_payload_bits.size()) * total_time_sec * 1.0e9);
      total_data_rates_gbps[index] = total_payload_bits / (total_time_sec * 1.0e9);
    }

    // Compute fired-channel occupancy using every raw chip instance as the denominator.
    if (threshold_sums[index].total_chip_instances > 0) {
      mean_active_channels[index] = static_cast<double>(threshold_sums[index].total_active_channels) /
          static_cast<double>(threshold_sums[index].total_chip_instances);
      chip_channel_rates_hz[index] = mean_active_channels[index] / kEventWindowSec;
    }
    
    // Store both requested percentiles for the three tail-rate products.
    for (std::size_t percentile = 0; percentile < kPercentiles.size(); ++percentile) {
      channel_hit_values[percentile][index] = channel_rates[percentile];
      chip_hit_values[percentile][index] = chip_rates[percentile];
      chip_data_values[percentile][index] = chip_data_rates[percentile];
    }
  }

  // Write channel, chip, payload, occupancy, and detector-total threshold scans.
  draw_overlay(parent->mkdir("channel_hit_rate_vs_mip"),
               "c_channel_hit_rate_vs_mip",
               "LFHCAL channel rates vs MIP threshold;MIP threshold;rate [Hz]",
               "g_channel_hit_rate_vs_mip",
               channel_hit_values);
  draw_overlay(parent->mkdir("chip_hit_rate_vs_mip"),
               "c_chip_hit_rate_vs_mip",
               "LFHCAL chip rates vs MIP threshold;MIP threshold;rate [Hz]",
               "g_chip_hit_rate_vs_mip",
               chip_hit_values);
  auto* chip_data_directory = parent->mkdir("chip_data_rate_vs_mip");
  draw_overlay(chip_data_directory,
               "c_chip_data_rate_vs_mip",
               "LFHCAL chip data rates vs MIP threshold;MIP threshold;data rate [Gb/s]",
               "g_chip_data_rate_vs_mip",
               chip_data_values);
  write_mean_graph(chip_data_directory, mean_chip_data_values);
  write_chip_channel_graphs(parent->mkdir("chip_channel_rate_vs_mip"),
                            chip_channel_rates_hz,
                            mean_active_channels);
  write_total_rate_graph(parent->mkdir("total_data_rate_vs_mip"), total_data_rates_gbps);
}

}  // namespace rates::lfhcal_analysis::occupancy::mip
