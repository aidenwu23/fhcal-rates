/*

./build/data_rate_vs_mip -i data/bkg_apr -o lfhcal_plots/data_rates/data_rate_vs_mip.root

*/

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {
// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
constexpr const char* kHitCollection = "LFHCALHits";
constexpr double kEventWindowSec = 2e-6;
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;
constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
constexpr std::array<double, 2> kPercentiles = {0.95, 0.99};

struct ThresholdProducts {
  // chip_bits[chip]: summed payload bits accumulated across the full sample.
  std::unordered_map<rates::LFHCALChipID, double, rates::LFHCALChipIDHash> chip_bits;
};

struct EventChip {
  // channel_energy[channel]: summed channel signal inside one chip for one event.
  std::unordered_map<rates::LFHCALChannelID, double, rates::LFHCALChannelIDHash> channel_energy;
};

// ----------------------------------------------------------------------------------
// CLI
// ----------------------------------------------------------------------------------
void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT_DIR -o OUTPUT.root\n";
}

struct Args {
  std::string input_dir;
  std::string output_file;
};

Args parse_args(int argc, char* argv[]) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_dir = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      args.output_file = argv[++i];
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input_dir.empty() || args.output_file.empty()) {
    usage(argv[0]);
    std::exit(1);
  }

  return args;
}

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

std::string percentile_title(double percentile) {
  return std::to_string(static_cast<int>(percentile * 100.0 + 0.5)) +
         "th percentile of data rate vs MIP coefficient;MIP coefficient;data rate [Gb/s]";
}


void draw_mean_graph(TFile& output,
                     const std::array<double, kCoefficients.size()>& values) {
  auto* dir = output.mkdir("mean");
  dir->cd();

  TCanvas canvas("c_mean_data_rate_vs_mip",
                 "Mean data rate vs MIP coefficient;MIP coefficient;data rate [Gb/s]",
                 1000,
                 800);
  canvas.SetGrid();
  canvas.SetLogy();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), 1.0e-4, kCoefficients.back(), 10.0);
  frame->SetTitle("Mean data rate vs MIP coefficient;MIP coefficient;data rate [Gb/s]");
  frame->SetStats(false);

  TGraph graph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), values.data());
  graph.SetName("g_mean_data_rate_vs_mip");
  graph.SetLineWidth(2);
  graph.SetMarkerStyle(20);
  graph.Draw("LP SAME");

  canvas.Write();
  graph.Write();
  output.cd();
}

void draw_percentile_graph(TFile& output,
                           double percentile,
                           const std::array<double, kCoefficients.size()>& values) {
  // Keep p95 and p99 products in separate top-level folders.
  auto* dir = output.mkdir(percentile_tag(percentile).c_str());
  dir->cd();

  TCanvas canvas(("c_" + percentile_tag(percentile) + "_data_rate_vs_mip").c_str(),
                 percentile_title(percentile).c_str(),
                 1000,
                 800);
  canvas.SetGrid();
  canvas.SetLogy();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), 1.0e-4, kCoefficients.back(), 10.0);
  frame->SetTitle(percentile_title(percentile).c_str());
  frame->SetStats(false);

  TGraph graph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), values.data());
  graph.SetName(("g_" + percentile_tag(percentile) + "_data_rate_vs_mip").c_str());
  graph.SetLineWidth(2);
  graph.SetMarkerStyle(20);
  graph.Draw("LP SAME");

  canvas.Write();
  graph.Write();
  output.cd();
}

void draw_overlay(TFile& output,
                  const std::array<double, kCoefficients.size()>& mean_values,
                  const std::array<std::array<double, kCoefficients.size()>, kPercentiles.size()>& percentile_values) {
  output.cd();

  TCanvas canvas("c_data_rate_vs_mip_overlay",
                 "Data rate vs MIP coefficient;MIP coefficient;data rate [Gb/s]",
                 1000,
                 800);
  canvas.SetGrid();
  canvas.SetLogy();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), 1.0e-4, kCoefficients.back(), 10.0);
  frame->SetTitle("Data rate vs MIP coefficient;MIP coefficient;data rate [Gb/s]");
  frame->SetStats(false);

  TLegend legend(0.58, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);

  TGraph g_mean(static_cast<int>(kCoefficients.size()), kCoefficients.data(), mean_values.data());
  g_mean.SetName("g_mean_data_rate_vs_mip_overlay");
  g_mean.SetLineWidth(2);
  g_mean.SetLineColor(kBlack);
  g_mean.SetMarkerColor(kBlack);
  g_mean.SetMarkerStyle(20);
  g_mean.Draw("LP SAME");
  legend.AddEntry(&g_mean, "mean", "lp");

  TGraph g95(static_cast<int>(kCoefficients.size()), kCoefficients.data(), percentile_values[0].data());
  g95.SetName("g_p95_data_rate_vs_mip_overlay");
  g95.SetLineWidth(2);
  g95.SetLineColor(kBlue + 1);
  g95.SetMarkerColor(kBlue + 1);
  g95.SetMarkerStyle(21);
  g95.Draw("LP SAME");
  legend.AddEntry(&g95, "95th percentile", "lp");

  TGraph g99(static_cast<int>(kCoefficients.size()), kCoefficients.data(), percentile_values[1].data());
  g99.SetName("g_p99_data_rate_vs_mip_overlay");
  g99.SetLineWidth(2);
  g99.SetLineColor(kRed + 1);
  g99.SetMarkerColor(kRed + 1);
  g99.SetMarkerStyle(22);
  g99.Draw("LP SAME");
  legend.AddEntry(&g99, "99th percentile", "lp");

  legend.Draw();
  canvas.Write();
  g_mean.Write();
  g95.Write();
  g99.Write();
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Grab args and files.
  const auto args = parse_args(argc, argv);
  const auto files = rates::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // products[threshold]: persistent bit sums for one MIP choice.
  std::array<ThresholdProducts, kCoefficients.size()> products{};
  const rates::LFHCALCellIDDecoder decoder;
  std::uint64_t n_events = 0;
  rates::FileProgress progress(files.size(), std::cerr);

  // Loop all files.
  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Loop events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!rates::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      // Channel sums are grouped under their parent chip.
      std::unordered_map<rates::LFHCALChipID, EventChip, rates::LFHCALChipIDHash> event_chips;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop all hits in this event.
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        // Increment the corresponding channel's energy for the corresponding readout chip.
        const auto channel = decoder.channel(cell_id);
        auto& event_chip = event_chips[decoder.decode_chip(cell_id)];
        event_chip.channel_energy[channel] += hit.getEnergy();
      }

      // After processing all hits into corresponding channels and chips, loop over all chips.
      for (const auto& [chip, event_chip] : event_chips) {
        std::array<int, kCoefficients.size()> active_counts{};

        // For each chip, loop through all channels.
        for (const auto& [channel, energy_gev] : event_chip.channel_energy) {
          const double mip_gev = rates::mip_energy_gev(channel.rlayerz);

          // For each channel, loop over all thresholds.
          for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

            // Skip if channel energy doesn't pass threshold.
            if (energy_gev <= kCoefficients[threshold_index] * mip_gev) continue;
            ++active_counts[threshold_index];
          }
        }

        // After summing all fired channels in this chip, convert it into data rate.
        // After summing all fired channels in this chip for this event, convert it into data rate.
        for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
          // Data size = (128 overhead bits + 32 bits per fired channel) * 4 samples
          if ((active_counts[threshold_index]) <= 0) continue;
          const double event_bits =
              (kOverheadBits + kBitsPerHit * static_cast<double>(active_counts[threshold_index])) * kSamplesPerEvent;
          products[threshold_index].chip_bits[chip] += event_bits;
        }
      }
    }
  }
  std::cerr << "\n";

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  // Convert accumulated bits into data rates, then extract mean and percentiles.
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<double, kCoefficients.size()> mean_values{};
  std::array<std::array<double, kCoefficients.size()>, kPercentiles.size()> percentile_values{};

  // Loop over all thresholds' accumulated products across the full sample.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    std::vector<double> rates_gbps;

    // Per threshold, loop through all accumulated chips within.
    for (const auto& [chip, total_bits] : products[threshold_index].chip_bits) {
      (void)chip;
      rates_gbps.push_back((total_bits / total_time_sec) / 1.0e9);
    }
    if (rates_gbps.empty()) continue;

    double sum = 0.0;
    for (const double rate_gbps : rates_gbps) {
      sum += rate_gbps;
    }
    mean_values[threshold_index] = sum / static_cast<double>(rates_gbps.size());

    // Sort so percentile lookup becomes a simple indexed read.
    std::sort(rates_gbps.begin(), rates_gbps.end());
    for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
      percentile_values[percentile_slot][threshold_index] =
          rates_gbps[percentile_index(rates_gbps.size(), kPercentiles[percentile_slot])];
    }
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  draw_mean_graph(output, mean_values);
  for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
    draw_percentile_graph(output, kPercentiles[percentile_slot], percentile_values[percentile_slot]);
  }
  draw_overlay(output, mean_values, percentile_values);

  output.Close();
  return 0;
}
