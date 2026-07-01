/*

./build/fired_channels -i data/bkg_apr -o plots/chip_occupancy/fired_channels.root

*/

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>

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
constexpr double kCoefficient = 0.5;

struct EventChip {
  // channel_energy[channel]: summed channel signal inside this chip for one event.
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
void draw_and_write(TDirectory* dir, TH1* hist, const char* canvas_name, bool logy = false) {
  dir->cd();
  TCanvas canvas(canvas_name, hist->GetTitle(), 1000, 800);
  if (logy) canvas.SetLogy();
  hist->SetStats(false);
  hist->Draw("hist");
  canvas.Write();
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse inputs and make sure there is data to process.
  const auto args = parse_args(argc, argv);
  const auto files = rates::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Histogram of how many channels fire in one chip per event.
  TH1D h_active_channels(
      "h_active_channels",
      "Fired channels per event per chip, 0.5 MIP threshold;fired channels/event/chip;Chips",
      28,
      0.0,
      28.0);

  const rates::LFHCALCellIDDecoder decoder;
  std::uint64_t n_events = 0;
  rates::FileProgress progress(files.size(), std::cerr);

  // Loop over all files.
  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Per file, loop over all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!rates::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      // For this event, group all channel sums under their parent chip.
      std::unordered_map<rates::LFHCALChipID, EventChip, rates::LFHCALChipIDHash> event_chips;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop hits.
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        // Increment the corresponding channel's energy for the corresponding readout chip.
        const auto channel = decoder.channel(cell_id);
        auto& event_chip = event_chips[decoder.decode_chip(cell_id)];

        event_chip.channel_energy[channel] += hit.getEnergy();
      }

      // active_counts[chip]: how many distinct channels passed threshold in this chip this event.
      std::unordered_map<rates::LFHCALChipID, int, rates::LFHCALChipIDHash> active_counts;

      // Loop over all chips.
      for (const auto& [chip, event_chip] : event_chips) {
        int active_count = 0;

        // Per chip, loop over all stored channels.
        for (const auto& [channel, energy_gev] : event_chip.channel_energy) {

          // Increment if that channel's energy exceeds the threshold.
          if (energy_gev > kCoefficient * rates::mip_energy_gev(channel.rlayerz)) ++active_count;
        }

        active_counts[chip] = active_count;
      }

      // Fill one histogram entry per chip seen in this event.
      for (const auto& [chip, count] : active_counts) {
        (void)chip;
        h_active_channels.Fill(count);
      }
    }
  }

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  // Write raw histograms and matching canvases.
  auto* hist_dir = output.mkdir("hists");
  hist_dir->cd();
  h_active_channels.Write();

  draw_and_write(&output, &h_active_channels, "c_active_channels", true);

  std::cout << "\n";
  output.Close();
  return 0;
}
