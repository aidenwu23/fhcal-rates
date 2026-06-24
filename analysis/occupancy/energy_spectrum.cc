/*

./build/energy_spectrum -i data/bkg_apr -o plots/occupancy/energy_spectrum.root

*/

#include <TCanvas.h>
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

// ----------------------------------------------------------------------------------
// CLI handling.
// ----------------------------------------------------------------------------------
struct Args {
  std::string input_dir;
  std::string output_file;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT_DIR -o OUTPUT.root\n";
}

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
void draw_and_write(TDirectory* dir, TH1* hist, const char* canvas_name, bool logx = false, bool logy = false) {
  dir->cd();
  TCanvas canvas(canvas_name, hist->GetTitle(), 1000, 800);
  if (logx) canvas.SetLogx();
  if (logy) canvas.SetLogy();
  hist->SetStats(false);
  hist->Draw("hist");
  canvas.Write();
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  const auto args = parse_args(argc, argv);
  const auto files = br::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  const auto energy_edges = br::log_edges(240, 1e-6, 10.0);
  std::array<TH1D*, kNReadoutLayers> h_channel_energy{};

  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    const std::string energy_name = "h_channel_energy_layer" + std::to_string(layer);
    const std::string energy_title =
        "LFHCAL channel energy per event, readout layer " + std::to_string(layer) +
        ";E_{dep}^{channel} [GeV];Channels";
    h_channel_energy[layer] = new TH1D(energy_name.c_str(), energy_title.c_str(), 240, energy_edges.data());
  }

  const br::LFHCALCellIDDecoder decoder;
  std::uint64_t n_events = 0;
  br::FileProgress progress(files.size(), std::cerr);

  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Histogrammed per-event? dunno if this works.
    // How to categorize a hit as being relevant to a channel other than being in the same event?
    // Two hits can come from the same channel but be very far away in time.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!br::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      std::array<std::unordered_map<br::LFHCALChannelID, double, br::LFHCALChannelIDHash>, kNReadoutLayers> channel_energy_by_layer;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop hits.
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        const auto channel = decoder.channel(cell_id);
        if (channel.rlayerz < 0 || channel.rlayerz >= kNReadoutLayers) continue;

        // Increment the corresponding channel's energy in the corresponding readout layer.
        channel_energy_by_layer[channel.rlayerz][channel] += hit.getEnergy();
      }

      // After accumualting all energy in each channel for this event, fill the corresponding hists.
      for (int layer = 0; layer < kNReadoutLayers; ++layer) {

        // Loop all filled channels.
        for (const auto& [channel, energy] : channel_energy_by_layer[layer]) {
          (void)channel;
          h_channel_energy[layer]->Fill(energy); // Fill corresponding layer.
        }
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

  auto* hist_dir = output.mkdir("hists");
  hist_dir->cd();
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    h_channel_energy[layer]->Write();
  }

  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    draw_and_write(&output, h_channel_energy[layer], ("c_layer_" + std::to_string(layer)).c_str(), true, true);
  }

  std::cout << "\n";
  return 0;
}
