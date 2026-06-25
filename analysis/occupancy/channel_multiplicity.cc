/*

./build/channel_multiplicity -i data/bkg_apr -o plots/occupancy/channel_multiplicity.root

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
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kTruthHitCollection = "LFHCALHits";
constexpr int kNLayers = 7;
using ThresholdsByLayer = std::array<double, kNLayers>;
constexpr double MIP_1 = 3.5e-3;
constexpr double MIP_2 = 7.0e-3;
constexpr double kCoefficient = 0.3;
const ThresholdsByLayer kThresholdsGeV = {
    kCoefficient * MIP_1,
    kCoefficient * MIP_1,
    kCoefficient * MIP_2,
    kCoefficient * MIP_2,
    kCoefficient * MIP_2,
    kCoefficient * MIP_2,
    kCoefficient * MIP_2};

// ----------------------------------------------------------------------------------
// CLI and per-event bookkeeping.
// ----------------------------------------------------------------------------------
struct Args {
  std::string input_dir;
  std::string output_file;
};

using EventCounts = std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>;

struct EventChannel {
  double energy_gev = 0.0;
  int hit_count = 0;
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
// Plot helpers.
// ----------------------------------------------------------------------------------
void draw_and_write(TDirectory* dir, TH1* hist, const char* canvas_name, bool logy = false) {
  dir->cd();
  TCanvas canvas(canvas_name, hist->GetTitle(), 1000, 800);
  if (logy) canvas.SetLogy();
  hist->SetStats(false);
  hist->Draw("hist");
  canvas.Write();
}

std::vector<TH1D*> make_layer_hists() {
  std::vector<TH1D*> hists;
  hists.reserve(kNLayers);
  for (int layer = 0; layer < kNLayers; ++layer) {
    const std::string name = "h_channel_multiplicity_layer" + std::to_string(layer);
    const std::string title = "LFHCAL channel multiplicity per event, layer " + std::to_string(layer) + ";hits/channel/event;Channels";
    hists.push_back(new TH1D(name.c_str(), title.c_str(), 100, 0.0, 100.0));
  }
  return hists;
}

void write_layer_hists(TDirectory* parent, const std::vector<TH1D*>& hists) {
  auto* hist_dir = parent->mkdir("hists");
  hist_dir->cd();
  for (auto* hist : hists) hist->Write();
  for (int layer = 0; layer < kNLayers; ++layer) {
    draw_and_write(parent, hists[layer], (std::string("c_layer") + std::to_string(layer)).c_str(), true);
  }
}

bool process_truth_event(const podio::Frame& frame,
                         const br::LFHCALCellIDDecoder& decoder,
                         const ThresholdsByLayer& thresholds_geV,
                         const std::vector<TH1D*>& hists) {
  if (!br::has_collection(frame, kTruthHitCollection)) return false;

  std::vector<EventCounts> event_counts(kNLayers);
  // event_channels[layer][channel]: For this event, summed channel signal before applying threshold.
  std::vector<std::unordered_map<br::LFHCALChannelID, EventChannel, br::LFHCALChannelIDHash>> event_channels(kNLayers);

  // Use LFHCALHits.
  const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kTruthHitCollection);
  for (const auto& hit : hits) {
    const auto channel_id = decoder.channel(static_cast<std::uint64_t>(hit.getCellID()));
    const int layer = channel_id.rlayerz;
    if (layer < 0 || layer >= kNLayers) continue;

    auto& event_channel = event_channels[layer][channel_id];
    event_channel.energy_gev += hit.getEnergy();
    ++event_channel.hit_count;
  }

  for (int layer = 0; layer < kNLayers; ++layer) {
    for (const auto& [channel_id, event_channel] : event_channels[layer]) {
      // Apply summed channel threshold.
      if (event_channel.energy_gev <= thresholds_geV[layer]) continue;
      event_counts[layer][channel_id] = event_channel.hit_count;
    }
  }

  for (int layer = 0; layer < kNLayers; ++layer) {
    for (const auto& [channel_id, count] : event_counts[layer]) {
      (void)channel_id;
      hists[layer]->Fill(count); // Fill hist.
    }
  }

  return true;
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

  const br::LFHCALCellIDDecoder decoder;
  auto hists = make_layer_hists();

  std::uint64_t n_events = 0;
  br::FileProgress progress(files.size(), std::cerr);

  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (process_truth_event(frame, decoder, kThresholdsGeV, hists)) ++n_events;
    }
  }

  if (n_events == 0) {
    std::cerr << "No events with LFHCALHits found\n";
    return 1;
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  write_layer_hists(&output, hists);

  std::cout << "\n";
  output.Close();
  return 0;
}
