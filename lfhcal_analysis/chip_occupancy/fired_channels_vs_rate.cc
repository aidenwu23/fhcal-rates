/*

./build/fired_channels_vs_rate -i data/bkg_apr -o lfhcal_plots/chip_occupancy/fired_channels_vs_rate.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "utils.h"

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
constexpr int kNReadoutLayers = 7;
constexpr double kCoefficient = 0.5;
constexpr double kEventWindowSec = 2e-6;
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;
constexpr std::array<double, 7> kContourRatesGbps = {0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035};
constexpr int kContourColor = 17;

struct ChipStats {
  std::uint64_t chip_passes = 0;
  std::uint64_t total_active_channels = 0;
};

struct EventChip {
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
// One active chip event contributes two pieces of serialized output:
// 1. Frame overhead.
// 2. Payload from fired channels.
double active_event_bits(double mean_fired_channels) {
  return (kOverheadBits + kBitsPerHit * mean_fired_channels) * kSamplesPerEvent;
}

// Convert one scatter-plot point into total chip output rate.
double total_chip_output_rate_gbps(double rate_hz, double mean_fired_channels) {
  return rate_hz * active_event_bits(mean_fired_channels) / 1.0e9;
}

// For one contour line, hold the total chip output rate fixed and solve for y.
double contour_y(double rate_hz, double rate_gbps) {
  return ((rate_gbps * 1.0e9 / rate_hz) / kSamplesPerEvent - kOverheadBits) / kBitsPerHit;
}

void draw_graph(TFile& output,
                const std::vector<double>& chip_rates_hz,
                const std::vector<double>& mean_fired_channels) {
  output.cd();

  TCanvas canvas("c_fired_channels_vs_rate",
                 "Mean fired channels per active event vs rate for each chip;rate [Hz];mean fired channels per active event",
                 1000,
                 800);
  canvas.SetGrid();

  constexpr double x_min = 0.0;
  constexpr double x_max = 5.0e4;
  constexpr double y_min = 0.0;
  constexpr double y_max = 5.0;
  auto* frame = canvas.DrawFrame(x_min, y_min, x_max, y_max);
  frame->SetTitle("Mean fired channels per active event vs rate for each chip;rate [Hz];mean fired channels per active event");
  frame->SetStats(false);

  TGraph graph(static_cast<int>(chip_rates_hz.size()), chip_rates_hz.data(), mean_fired_channels.data());
  graph.SetName("g_fired_channels_vs_rate");
  graph.SetMarkerStyle(20);
  graph.Draw("P SAME");

  std::vector<TGraph> contour_graphs;
  contour_graphs.reserve(kContourRatesGbps.size());

  // Build one curve per target Gb/s value.
  for (std::size_t contour_index = 0; contour_index < kContourRatesGbps.size(); ++contour_index) {
    std::vector<double> x_values;
    std::vector<double> y_values;

    // Step across the visible x range, then solve for the matching y value.
    for (int step = 1; step <= 400; ++step) {
      const double rate_hz = x_min + (x_max - x_min) * static_cast<double>(step) / 400.0;
      const double y_value = contour_y(rate_hz, kContourRatesGbps[contour_index]);

      // Keep the point only when the solved y value falls inside the visible frame.
      if (y_value < y_min || y_value > y_max) continue;
      x_values.push_back(rate_hz);
      y_values.push_back(y_value);
    }
    if (x_values.empty()) continue;

    // Convert the sampled x/y points into one drawable ROOT curve.
    contour_graphs.emplace_back(static_cast<int>(x_values.size()), x_values.data(), y_values.data());
    auto& contour = contour_graphs.back();
    contour.SetName(("g_contour_" + std::to_string(contour_index)).c_str());
    contour.SetLineColor(kContourColor);
    contour.SetLineStyle(2);
    contour.SetLineWidth(2);
    contour.Draw("L SAME");
  }


  canvas.Write();
  graph.Write();
  for (auto& contour : contour_graphs) contour.Write();
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Handle args and input files.
  const auto args = parse_args(argc, argv);
  const auto files = rates::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  const rates::LFHCALCellIDDecoder decoder;
  std::uint64_t n_events = 0;
  std::unordered_map<rates::LFHCALChipID, ChipStats, rates::LFHCALChipIDHash> chip_stats;
  rates::FileProgress progress(files.size(), std::cerr);

  // Loop all files.
  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Per file, loop all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!rates::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      std::unordered_map<rates::LFHCALChipID, EventChip, rates::LFHCALChipIDHash> event_chips;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Per event, loop all hits.
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
        int active_count = 0;

        // Per chip, loop over all channels.
        for (const auto& [channel, energy_gev] : event_chip.channel_energy) {

          // Increment if channel exceeds threshold.
          if (energy_gev > kCoefficient * rates::mip_energy_gev(channel.rlayerz)) ++active_count;
        }
        if (active_count <= 0) continue;

        auto& stats = chip_stats[chip];
        stats.total_active_channels += static_cast<std::uint64_t>(active_count);
        ++stats.chip_passes;
      }
    }
  }
  std::cerr << "\n";

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::vector<double> chip_rates_hz;
  std::vector<double> mean_fired_channels;
  chip_rates_hz.reserve(chip_stats.size());
  mean_fired_channels.reserve(chip_stats.size());

  // After collecting data from all events, process the statistics pooled in each chip.
  for (const auto& [chip, stats] : chip_stats) {
    (void)chip;
    chip_rates_hz.push_back(static_cast<double>(stats.chip_passes) / total_time_sec);
    // For the y axis, count only events where the chip actually fires.
    mean_fired_channels.push_back(static_cast<double>(stats.total_active_channels) / static_cast<double>(stats.chip_passes));
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  draw_graph(output, chip_rates_hz, mean_fired_channels);

  output.Close();
  return 0;
}
