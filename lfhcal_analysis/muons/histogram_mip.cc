/*

./build/histogram_mip -i data/mu-_10GeV_lfhcal_100k.edm4hep.root -o lfhcal_plots/muons/histogram_mip.root

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

constexpr const char* kHitCollection = "LFHCALHits";
constexpr int kNReadoutLayers = 7;
constexpr double kMaxRadiusMm = 1000.0;
constexpr double kMaxRadiusMm2 = kMaxRadiusMm * kMaxRadiusMm;
constexpr double kMinContributionGeV = 0.0005;

// ----------------------------------------------------------------------------------
// CLI handling.
// ----------------------------------------------------------------------------------
struct Args {
  std::string input_file;
  std::string output_file;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT.root -o OUTPUT.root\n";
}

// Parse CLI arguments.
Args parse_args(int argc, char* argv[]) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_file = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      args.output_file = argv[++i];
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input_file.empty() || args.output_file.empty()) {
    usage(argv[0]);
    std::exit(1);
  }

  return args;
}

// Thing.
void draw_and_write(TDirectory* dir, TH1* hist, const std::string& canvas_name) {
  dir->cd();
  TCanvas canvas(canvas_name.c_str(), hist->GetTitle(), 1000, 800);
  canvas.SetLogy();
  hist->Draw("hist");
  rates::pad_axes(*hist);
  canvas.Write();
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and create the output directory if needed.
  const auto args = parse_args(argc, argv);
  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Read the single input ROOT file.
  podio::ROOTReader reader;
  reader.openFile(args.input_file);
  const std::size_t total_events = reader.getEntries("events");
  if (total_events == 0) {
    std::cerr << "No events found in " << args.input_file << "\n";
    return 1;
  }

  // All physical layers use the same tile, so one histogram collects every tile energy.
  auto* h_tile_energy = new TH1D(
      "h_tile_energy",
      "Muon tile energy per event;tile energy [GeV];counts",
      200,
      0.0,
      0.005);

  // Channel construction changes with readout layer, so keep one histogram per layer.
  std::array<TH1D*, kNReadoutLayers> h_channel_energy{};
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    h_channel_energy[layer] = new TH1D(
        ("h_channel_energy_layer" + std::to_string(layer)).c_str(),
        ("Muon channel energy per event, layer " + std::to_string(layer) + ";channel energy [GeV];counts").c_str(),
        200,
        0.0,
        0.05);
  }

  const rates::LFHCALCellIDDecoder decoder;

  // Loop events.
  for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
    auto data = reader.readEntry("events", event_index);
    if (!data) continue;

    podio::Frame frame(std::move(data));
    if (!rates::has_collection(frame, kHitCollection)) continue;

    // For this event, sum each physical tile and each layer-dependent readout channel.
    std::unordered_map<std::uint64_t, double> tile_energy;
    std::array<std::unordered_map<rates::LFHCALChannelID, double, rates::LFHCALChannelIDHash>, kNReadoutLayers> channel_energy_by_layer;
    const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

    // Loop hits and assign each accepted energy deposit to its tile and channel sums.
    for (const auto& hit : hits) {
      const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
      if (decoder.is_passive(cell_id)) continue;

      const auto channel = decoder.channel(cell_id);
      if (channel.rlayerz < 0 || channel.rlayerz >= kNReadoutLayers) continue;

      // Keep tiles and channels whose transverse radius is below 1000 mm.
      const auto position = decoder.position(cell_id);
      const double radius_mm2 = position.x_mm * position.x_mm + position.y_mm * position.y_mm;
      if (radius_mm2 > kMaxRadiusMm2) continue;

      double filtered_hit_energy_gev = 0.0;

      // Skip tiny contributions before adding this hit to the tile and channel sums.
      for (const auto& contribution : hit.getContributions()) {
        if (contribution.getEnergy() < kMinContributionGeV) continue;
        filtered_hit_energy_gev += contribution.getEnergy();
      }

      if (filtered_hit_energy_gev > 0.0) {
        tile_energy[cell_id] += filtered_hit_energy_gev;
        channel_energy_by_layer[channel.rlayerz][channel] += filtered_hit_energy_gev;
      }
    }

    // Once the full event has been summed, fill one entry per tile.
    for (const auto& [cell_id, energy_gev] : tile_energy) {
      (void)cell_id;
      h_tile_energy->Fill(energy_gev);
    }

    // Fill channels separately because their summed depth depends on the readout layer.
    for (int layer = 0; layer < kNReadoutLayers; ++layer) {
      for (const auto& [channel, energy_gev] : channel_energy_by_layer[layer]) {
        (void)channel;
        h_channel_energy[layer]->Fill(energy_gev);
      }
    }
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  auto* hist_dir = output.mkdir("hists");
  hist_dir->cd();
  h_tile_energy->Write();
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    h_channel_energy[layer]->Write();
  }

  draw_and_write(&output, h_tile_energy, "c_tile_energy");
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    draw_and_write(&output,
                   h_channel_energy[layer],
                   "c_channel_energy_layer" + std::to_string(layer));
  }

  output.Close();
  return 0;
}
