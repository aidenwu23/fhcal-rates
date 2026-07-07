/*

./build/lfhcal_rate -i data/bkg_apr -o lfhcal_plots/data_rates/lfhcal_rate.root

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
constexpr double kEventWindowSec = 2e-6;
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;
constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

struct ThresholdProducts {
  double total_bits = 0.0;
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
void draw_graph(TFile& output, const std::array<double, kCoefficients.size()>& rates_gbps) {
  output.cd();

  TCanvas canvas("c_lfhcal_rate",
                 "LFHCAL total data rate vs MIP coefficient;MIP coefficient;data rate [Gb/s]",
                 1000,
                 800);
  canvas.SetGrid();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), 0.0, kCoefficients.back(), 55.0);
  frame->SetTitle("LFHCAL total data rate vs MIP coefficient;MIP coefficient;data rate [Gb/s]");
  frame->SetStats(false);

  TGraph graph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), rates_gbps.data());
  graph.SetName("g_lfhcal_rate");
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

  // Inputs.
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

  // Loop all input files.
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

      // Per event, loop thru all hits.
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
        std::array<int, kCoefficients.size()> active_counts{};

        // For each chip, loop over all of it's channels.
        for (const auto& [channel, energy_gev] : event_chip.channel_energy) {
          const double mip_gev = rates::mip_energy_gev(channel.rlayerz);

          // For each channel, loop through all MIP thresholds.
          for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
            if (energy_gev <= kCoefficients[threshold_index] * mip_gev) continue;

            // Increment if channel passes this threshold.
            ++active_counts[threshold_index];
          }
        }

        // After computing which channels pass for which thresholds, compute the corresponding data rate
        // at each threshold and add it onto the total data rate for the whole lfhcal.
        // After summing all fired channels in this chip for this event, convert it into data rate.
        for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
          if ((active_counts[threshold_index]) <= 0) continue;
          const double event_bits =
              (kOverheadBits + kBitsPerHit * static_cast<double>(active_counts[threshold_index])) * kSamplesPerEvent;
          products[threshold_index].total_bits += event_bits;
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
  std::array<double, kCoefficients.size()> rates_gbps{};
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    rates_gbps[threshold_index] = (products[threshold_index].total_bits / total_time_sec) / 1.0e9;
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  draw_graph(output, rates_gbps);

  output.Close();
  return 0;
}
