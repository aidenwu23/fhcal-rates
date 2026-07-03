/*

./build/chip_channel_rate_vs_mip -i data/bkg_apr -o plots/chip_occupancy/chip_channel_rate_vs_mip.root

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

namespace fs = std::filesystem;

namespace {
// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
constexpr const char* kHitCollection = "LFHCALHits";
constexpr int kNReadoutLayers = 7;
constexpr double kEventWindowSec = 2e-6;
constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

struct ThresholdProducts {
  std::uint64_t total_active_channels = 0;
  std::uint64_t total_chip_instances = 0;
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
void draw_graph(TFile& output, const std::array<double, kCoefficients.size()>& rates_hz) {
  output.cd();

  TCanvas canvas("c_chip_channel_rate_vs_mip",
                 "Mean fired channel rate per chip vs MIP coefficient;MIP coefficient;fired channel rate [Hz/chip]",
                 1000,
                 800);
  canvas.SetGrid();
  canvas.SetLogy();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), 1.0e4, kCoefficients.back(), 1.0e7);
  frame->SetTitle("Mean fired channel rate per chip vs MIP coefficient;MIP coefficient;fired channel rate [Hz/chip]");
  frame->SetStats(false);

  TGraph graph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), rates_hz.data());
  graph.SetName("g_chip_channel_rate_vs_mip");
  graph.SetLineWidth(2);
  graph.SetMarkerStyle(20);
  graph.Draw("LP SAME");

  canvas.Write();
  graph.Write();
}

void draw_active_channels_graph(TFile& output, const std::array<double, kCoefficients.size()>& means) {
  output.cd();

  TCanvas canvas("c_active_channels_vs_mip",
                 "Mean fired channels per event per chip vs MIP coefficient;MIP coefficient;mean fired channels/event/chip",
                 1000,
                 800);
  canvas.SetGrid();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), 0.0, kCoefficients.back(), 128.0);
  frame->SetTitle("Mean fired channels per event per chip vs MIP coefficient;MIP coefficient;mean fired channels/event/chip");
  frame->SetStats(false);

  TGraph graph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), means.data());
  graph.SetName("g_active_channels_vs_mip");
  graph.SetLineWidth(2);
  graph.SetMarkerStyle(20);
  graph.Draw("LP SAME");

  canvas.Write();
  graph.Write();
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Grab args and input events.
  const auto args = parse_args(argc, argv);
  const auto files = rates::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

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
        (void)chip;

        // For each chip, loop over all thresholds.
        for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
          int active_count = 0;

          // Per threshold, loop over all channels in that chip.
          for (const auto& [channel, energy_gev] : event_chip.channel_energy) {
            
            // Increment if a channel is above a threshold.
            if (energy_gev > kCoefficients[threshold_index] * rates::mip_energy_gev(channel.rlayerz)) ++active_count;
          }
          products[threshold_index].total_active_channels += static_cast<std::uint64_t>(active_count);
          ++products[threshold_index].total_chip_instances;
        }
      }
    }
  }
  std::cerr << "\n";

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  std::array<double, kCoefficients.size()> rates_hz{};
  std::array<double, kCoefficients.size()> means{};
  // Loop over all thresholds' accumulated products across the full sample.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    const auto& product = products[threshold_index];
    if (product.total_chip_instances == 0) continue;
    means[threshold_index] =
        static_cast<double>(product.total_active_channels) / static_cast<double>(product.total_chip_instances);
    rates_hz[threshold_index] =
        static_cast<double>(product.total_active_channels) /
        (static_cast<double>(product.total_chip_instances) * kEventWindowSec);
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  draw_graph(output, rates_hz);
  draw_active_channels_graph(output, means);

  output.Close();
  return 0;
}
