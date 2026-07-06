/*

./build/data_rate_vs_radius -i data/bkg_apr -o plots/data_rates/data_rate_vs_radius.root

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
#include <cmath>
#include <cstdint>
#include <cstdio>
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
constexpr std::array<double, 3> kCoefficients = {0.1, 0.5, 1.5};
constexpr std::array<double, 3> kPercentiles = {0.0, 0.95, 0.99};
constexpr int kRadiusBins = 30;
constexpr double kRadiusMaxMm = 3000.0;
const std::array<int, kPercentiles.size()> kColors = {kBlack, kBlue + 1, kRed + 1};
const std::array<int, kPercentiles.size()> kMarkers = {20, 21, 22};

struct ThresholdProducts {
  // chip_bits[chip]: summed payload bits accumulated across the full sample.
  std::unordered_map<rates::LFHCALChipID, double, rates::LFHCALChipIDHash> chip_bits;
};

struct EventChip {
  // channel_energy[channel]: summed channel signal inside one chip for one event.
  std::unordered_map<rates::LFHCALChannelID, double, rates::LFHCALChannelIDHash> channel_energy;
};

struct RadiusBinStats {
  std::vector<double> rates_gbps;
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
std::string threshold_label(double threshold_mip) {
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "%.2f MIP", threshold_mip);
  return buffer;
}

std::string threshold_tag(double threshold_mip) {
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "thr%03dMIP", static_cast<int>(threshold_mip * 100.0 + 0.5));
  return buffer;
}

std::string percentile_name(double percentile) {
  if (percentile == 0.0) return "mean";
  return "p" + std::to_string(static_cast<int>(percentile * 100.0 + 0.5));
}

std::string percentile_label(double percentile) {
  if (percentile == 0.0) return "mean";
  return std::to_string(static_cast<int>(percentile * 100.0 + 0.5)) + "th percentile";
}

// Turn a percentile request into the corresponding sorted-data index.
std::size_t percentile_index(std::size_t n_values, double percentile) {
  if (n_values == 0) return 0;
  return static_cast<std::size_t>(percentile * static_cast<double>(n_values - 1));
}



std::size_t radius_bin(double radius_mm_value) {
  if (radius_mm_value >= kRadiusMaxMm) return kRadiusBins - 1;
  const double fraction = radius_mm_value / kRadiusMaxMm;
  return std::min<std::size_t>(kRadiusBins - 1, static_cast<std::size_t>(fraction * kRadiusBins));
}

void draw_overlay(TFile& output,
                  double threshold_mip,
                  const std::array<std::array<double, kRadiusBins>, kPercentiles.size()>& values,
                  const std::array<double, kRadiusBins>& radius_centers_mm) {
  // Keep one radius scan per threshold in its own top-level folder.
  auto* dir = output.mkdir(threshold_tag(threshold_mip).c_str());
  dir->cd();

  TCanvas canvas(("c_" + threshold_tag(threshold_mip) + "_data_rate_vs_radius").c_str(),
                 ("Chip data rate vs radius, threshold " + threshold_label(threshold_mip) + ";radius [mm];chip data rate [Gb/s]").c_str(),
                 1000,
                 800);
  canvas.SetGrid();
  canvas.SetLogy();

  auto* frame = canvas.DrawFrame(350.0, 1.0e-6, 2850.0, 4.0e-2);
  frame->SetTitle(("Chip data rate vs radius, threshold " + threshold_label(threshold_mip) + ";radius [mm];chip data rate [Gb/s]").c_str());
  frame->SetStats(false);

  TLegend legend(0.62, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);

  std::array<TGraph, kPercentiles.size()> graphs = {
      TGraph(kRadiusBins, radius_centers_mm.data(), values[0].data()),
      TGraph(kRadiusBins, radius_centers_mm.data(), values[1].data()),
      TGraph(kRadiusBins, radius_centers_mm.data(), values[2].data())};

  for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
    graphs[percentile_slot].SetName(("g_" + percentile_name(kPercentiles[percentile_slot]) + "_" + threshold_tag(threshold_mip) + "_data_rate_vs_radius").c_str());
    graphs[percentile_slot].SetLineColor(kColors[percentile_slot]);
    graphs[percentile_slot].SetMarkerColor(kColors[percentile_slot]);
    graphs[percentile_slot].SetMarkerStyle(kMarkers[percentile_slot]);
    graphs[percentile_slot].SetLineWidth(2);
    graphs[percentile_slot].Draw("LP SAME");
    legend.AddEntry(&graphs[percentile_slot], percentile_label(kPercentiles[percentile_slot]).c_str(), "lp");
  }

  legend.Draw();
  canvas.Write();
  for (auto& graph : graphs) graph.Write();
  output.cd();
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

  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<double, kRadiusBins> radius_centers_mm{};
  const double radius_bin_width_mm = kRadiusMaxMm / static_cast<double>(kRadiusBins);
  for (int bin = 0; bin < kRadiusBins; ++bin) {
    radius_centers_mm[bin] = (static_cast<double>(bin) + 0.5) * radius_bin_width_mm;
  }

  std::array<std::array<std::array<double, kRadiusBins>, kPercentiles.size()>, kCoefficients.size()> plot_values{};

  // Loop over all thresholds.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    std::array<RadiusBinStats, kRadiusBins> bins{};

    // For each threshold, loop over all chips.
    for (const auto& [chip, total_bits] : products[threshold_index].chip_bits) {

      // Compute stats for this chip at this threshold and put it in the corresponding radius bin.
      const double rate_gbps = (total_bits / total_time_sec) / 1.0e9;
      const auto chip_position = decoder.position(chip);
      const double chip_radius_mm = std::hypot(chip_position.x_mm, chip_position.y_mm);
      bins[radius_bin(chip_radius_mm)].rates_gbps.push_back(rate_gbps);
    }

    // For this threshold, loop through all radius bins.
    for (int bin = 0; bin < kRadiusBins; ++bin) {
      auto& rates_gbps = bins[bin].rates_gbps;
      if (rates_gbps.empty()) continue;

      double sum = 0.0;
      for (const double rate_gbps : rates_gbps) sum += rate_gbps;
      plot_values[threshold_index][0][bin] = sum / static_cast<double>(rates_gbps.size());

      // Sort so percentile lookup becomes a simple indexed read.
      std::sort(rates_gbps.begin(), rates_gbps.end());
      plot_values[threshold_index][1][bin] = rates_gbps[percentile_index(rates_gbps.size(), 0.95)];
      plot_values[threshold_index][2][bin] = rates_gbps[percentile_index(rates_gbps.size(), 0.99)];
    }
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    draw_overlay(output, kCoefficients[threshold_index], plot_values[threshold_index], radius_centers_mm);
  }

  output.Close();
  return 0;
}
