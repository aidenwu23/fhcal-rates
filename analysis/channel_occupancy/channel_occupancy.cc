/*

./build/channel_occupancy -i data/bkg_apr -o plots/channel_occupancy/channel_occupancy.root

*/

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include "shared.h"
#include "truth.h"
#include "utils.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

namespace {

// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
constexpr double MIP_1 = 3.5e-3;
constexpr double MIP_2 = 7.0e-3;
constexpr double kCoefficient = 0.5;
const rates::channel_occupancy::ThresholdsByLayer kThresholdsGeV = {
    kCoefficient * MIP_1,
    kCoefficient * MIP_1,
    kCoefficient * MIP_2,
    kCoefficient * MIP_2,
    kCoefficient * MIP_2,
    kCoefficient * MIP_2,
    kCoefficient * MIP_2};

struct Args {
  std::string input_dir;
  std::string output_file;
};

// ----------------------------------------------------------------------------------
// CLI
// ----------------------------------------------------------------------------------
void usage(const char* argv0);

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

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT_DIR -o OUTPUT.root\n";
}

// ----------------------------------------------------------------------------------
// Plot helper(s)
// ----------------------------------------------------------------------------------
void write_layer_directory(TDirectory* parent,
                           const rates::channel_occupancy::LayerAccum& layer_accum,
                           int layer,
                           std::uint64_t n_events,
                           const std::optional<rates::channel_occupancy::AxisEdges2D>& reference_edges = std::nullopt) {
  const std::string dir_name = layer == rates::channel_occupancy::kAllLayersIndex ? "sum_layers" : std::string("layer") + std::to_string(layer);
  auto* dir = parent->mkdir(dir_name.c_str());
  dir->cd();

  if (layer_accum.channels.empty()) {
    auto* hist_dir = dir->mkdir("hists");
    hist_dir->cd();
    layer_accum.h_hits_evt->Write();
    return;
  }

  double rmax = 0.0;
  for (const auto& [channel_id, stats] : layer_accum.channels) {
    (void)channel_id;
    rmax = std::max(rmax, stats.r());
  }
  const rates::channel_occupancy::AxisEdges2D axis_edges =
      reference_edges.has_value() ? *reference_edges : rates::channel_occupancy::make_layer_axis_edges(layer_accum);
  const auto& x_edges = axis_edges.x_edges;
  const auto& y_edges = axis_edges.y_edges;

  auto* h_avg = new TH2D(
      "h_avg",
      (layer == rates::channel_occupancy::kAllLayersIndex ? std::string("Summed layers;x [mm];y [mm];avg hits/event/channel")
        : std::string("Layer ") + std::to_string(layer) + ";x [mm];y [mm];avg hits/event/channel").c_str(),
      static_cast<int>(x_edges.size()) - 1, x_edges.data(),
      static_cast<int>(y_edges.size()) - 1, y_edges.data());

  auto* h_rate = new TH2D(
      "h_rate",
      (layer == rates::channel_occupancy::kAllLayersIndex ? std::string("Summed layers;x [mm];y [mm];rate [Hz/channel]")
        : std::string("Layer ") + std::to_string(layer) + ";x [mm];y [mm];rate [Hz/channel]").c_str(),
      static_cast<int>(x_edges.size()) - 1, x_edges.data(),
      static_cast<int>(y_edges.size()) - 1, y_edges.data());

  auto* h_avg_r = new TH1D(
      "h_avg_r",
      (layer == rates::channel_occupancy::kAllLayersIndex ? std::string("Summed layers;R [mm];avg hits/event")
        : std::string("Layer ") + std::to_string(layer) + ";R [mm];avg hits/event").c_str(),
      100, 0.0, std::max(1.0, 1.05 * rmax));

  auto* h_nchan_r = new TH1D("h_nchan_r", "", 100, 0.0, std::max(1.0, 1.05 * rmax));

  // Fill hists.
  for (const auto& [channel_id, stats] : layer_accum.channels) {
    (void)channel_id;
    const double avg = static_cast<double>(stats.total_hits) / static_cast<double>(n_events);
    const double rate = static_cast<double>(stats.total_hits) / (static_cast<double>(n_events) * rates::channel_occupancy::kEventWindowSec);
    h_avg->Fill(stats.x_mm, stats.y_mm, avg);
    h_rate->Fill(stats.x_mm, stats.y_mm, rate);
    h_avg_r->Fill(stats.r(), avg);
    h_nchan_r->Fill(stats.r(), 1.0);
  }

  for (int bin = 1; bin <= h_avg_r->GetNbinsX(); ++bin) {
    const double nch = h_nchan_r->GetBinContent(bin);
    if (nch > 0.0) h_avg_r->SetBinContent(bin, h_avg_r->GetBinContent(bin) / nch);
  }

  auto* hist_dir = dir->mkdir("hists");
  hist_dir->cd();
  h_avg->Write();
  h_rate->Write();
  h_avg_r->Write();
  h_nchan_r->Write();
  layer_accum.h_hits_evt->Write();

  rates::channel_occupancy::draw_and_write(dir, h_avg, "c_avg", true, false);
  rates::channel_occupancy::draw_and_write(dir, h_rate, "c_rate", true, false);
  rates::channel_occupancy::draw_and_write(dir, h_avg_r, "c_avg_r", false, false);
  rates::channel_occupancy::draw_and_write(dir, layer_accum.h_hits_evt, "c_hits_evt", false, true);
}

// Truth has access to contributions, so write one for each origin type in addition to one for every layer.
void write_truth(TFile& output,
                 const std::vector<rates::channel_occupancy::TruthOccupancyGroup>& truth_groups,
                 std::uint64_t n_events) {

  // Use DIS for bin edges since DIS typically has lots of stats.
  std::vector<std::optional<rates::channel_occupancy::AxisEdges2D>> dis_reference_edges(rates::channel_occupancy::kNLayers + 1);
  const auto dis_it = std::find_if(
      truth_groups.begin(),
      truth_groups.end(),
      [](const rates::channel_occupancy::TruthOccupancyGroup& group) { return group.label == "DIS"; });
  if (dis_it != truth_groups.end()) {
    for (int layer = 0; layer <= rates::channel_occupancy::kNLayers; ++layer) {
      if (!dis_it->layers[layer].channels.empty()) {
        dis_reference_edges[layer] = rates::channel_occupancy::make_layer_axis_edges(dis_it->layers[layer]);
      }
    }
  }

  for (const auto& truth_group : truth_groups) {
    auto* parent = output.mkdir(truth_group.label.c_str());
    for (int layer = 0; layer <= rates::channel_occupancy::kNLayers; ++layer) {
      write_layer_directory(parent, truth_group.layers[layer], layer, n_events, dis_reference_edges[layer]);
    }
  }
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and initialize pointers.
  const auto args = parse_args(argc, argv);
  const auto& input_dir = args.input_dir;
  const auto& output_file = args.output_file;
  const auto& thresholds_geV = kThresholdsGeV;

  // Create a decoder for cell IDs.
  const rates::LFHCALCellIDDecoder decoder;
  // Find input files.
  const auto files = rates::find_root_files(input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << input_dir << "\n";
    return 1;
  }

  // Ensure output directory.
  fs::path output_path = output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // NOTE: reco unavailable for now due to some timing issues with background.
  std::vector<rates::channel_occupancy::TruthOccupancyGroup> truth_groups;
  rates::channel_occupancy::init_truth_groups(truth_groups);

  std::uint64_t n_truth_events = 0;

  rates::FileProgress progress(files.size(), std::cerr);

  // Loop over all files in the input directory.
  for (const auto& path : files) {
    progress.tick();

    podio::ROOTReader reader;
    reader.openFile(path.string());

    const std::size_t total_events = reader.getEntries("events");

    // Per file, loop over all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {

      // Read event.
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));

      if (rates::channel_occupancy::process_truth_event(frame, decoder, thresholds_geV, truth_groups)) ++n_truth_events;
    }
  }

  if (n_truth_events == 0) {
    std::cerr << "No events with LFHCALHits found\n";
    return 1;
  }

  TFile output(output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << output_file << "\n";
    return 1;
  }

  if (n_truth_events > 0) write_truth(output, truth_groups, n_truth_events);

  std::cout << "\n";
  output.Close();
  return 0;
}
