/*

./build/occupancy -i data/reco_bkg_apr -o plots/occupancy/occupancy.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include "decode_channel.h"
#include "reco.h"
#include "smooth_hists.h"
#include "truth.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

namespace {

// Constants
constexpr double kDefaultThresholdGeV = 5e-4;           // Default threshold is 0.5 MeV.
constexpr double kEventWindowSec = 2e-6;                // 2 microseconds.
constexpr double kCellSizeMM = 51.0;                    // Display bin size in x and y.
constexpr double kXYExtentMM = 2700;                    // +- range for x and y.

struct Args {
  std::string input_dir;
  std::string output_file;
  double threshold_geV = kDefaultThresholdGeV;
};

void usage(const char* argv0);

// Parse CLI arguments.
Args parse_args(int argc, char* argv[]) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_dir = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      args.output_file = argv[++i];
    } else if ((arg == "-t" || arg == "--threshold") && i + 1 < argc) {
      args.threshold_geV = std::stod(argv[++i]);
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
  std::cerr << "Usage: " << argv0 << " -i INPUT_DIR -o OUTPUT.root [-t THRESHOLD_GEV]\n";
}

// Draw a (styled) histogram on a canvas and write it to the given directory.
void draw_and_write(TDirectory* canvas_dir, TH1* hist, const char* canvas_name, bool logz = false, bool logy = false) {
  canvas_dir->cd();
  TCanvas canvas(canvas_name, hist->GetTitle(), 1000, 800);
  if (logz) canvas.SetLogz();
  if (logy) canvas.SetLogy();
  hist->SetStats(false);
  const bool is2d = hist->InheritsFrom(TH2::Class());
  hist->Draw(is2d ? "colz" : "hist");
  canvas.Write();
}

void write_layer_directory(TDirectory* parent,
                           const br::occupancy::LayerAccum& layer_accum,
                           int layer,
                           std::uint64_t n_events) {
  const std::string dir_name = layer == br::occupancy::kAllLayersIndex ? "sum_layers" : std::string("layer") + std::to_string(layer);
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
  const double xmin = -kXYExtentMM;
  const double xmax = kXYExtentMM;
  const double ymin = -kXYExtentMM;
  const double ymax = kXYExtentMM;
  const int xbins = static_cast<int>(std::ceil((xmax - xmin) / kCellSizeMM));
  const int ybins = static_cast<int>(std::ceil((ymax - ymin) / kCellSizeMM));

  auto* h_avg = new TH2D(
      "h_avg",
      (layer == br::occupancy::kAllLayersIndex ? std::string("Summed layers;x [mm];y [mm];avg hits/event/channel")
        : std::string("Layer ") + std::to_string(layer) + ";x [mm];y [mm];avg hits/event/channel").c_str(),
      xbins, xmin, xmax,
      ybins, ymin, ymax);

  auto* h_max = new TH2D(
      "h_max",
      (layer == br::occupancy::kAllLayersIndex ? std::string("Summed layers;x [mm];y [mm];max hits/event/channel")
        : std::string("Layer ") + std::to_string(layer) + ";x [mm];y [mm];max hits/event/channel").c_str(),
      xbins, xmin, xmax,
      ybins, ymin, ymax);

  auto* h_rate = new TH2D(
      "h_rate",
      (layer == br::occupancy::kAllLayersIndex ? std::string("Summed layers;x [mm];y [mm];rate [Hz/channel]")
        : std::string("Layer ") + std::to_string(layer) + ";x [mm];y [mm];rate [Hz/channel]").c_str(),
      xbins, xmin, xmax,
      ybins, ymin, ymax);

  auto* h_avg_r = new TH1D(
      "h_avg_r",
      (layer == br::occupancy::kAllLayersIndex ? std::string("Summed layers;R [mm];avg hits/event")
        : std::string("Layer ") + std::to_string(layer) + ";R [mm];avg hits/event").c_str(),
      100, 0.0, std::max(1.0, 1.05 * rmax));

  auto* h_nchan_r = new TH1D("h_nchan_r", "", 100, 0.0, std::max(1.0, 1.05 * rmax));

  for (const auto& [channel_id, stats] : layer_accum.channels) {
    (void)channel_id;
    const double avg = static_cast<double>(stats.total_hits) / static_cast<double>(n_events);
    const double rate = static_cast<double>(stats.total_hits) / (static_cast<double>(n_events) * kEventWindowSec);
    h_avg->Fill(stats.x(), stats.y(), avg);
    h_max->Fill(stats.x(), stats.y(), stats.max_hits_event);
    h_rate->Fill(stats.x(), stats.y(), rate);
    h_avg_r->Fill(stats.r(), avg);
    h_nchan_r->Fill(stats.r(), 1.0);
  }

  for (int bin = 1; bin <= h_avg_r->GetNbinsX(); ++bin) {
    const double nch = h_nchan_r->GetBinContent(bin);
    if (nch > 0.0) h_avg_r->SetBinContent(bin, h_avg_r->GetBinContent(bin) / nch);
  }

  br::smooth_hist_vertical(h_avg);
  br::smooth_hist_vertical(h_max);
  br::smooth_hist_vertical(h_rate);
  br::smooth_hist_neighbhors(h_avg);
  br::smooth_hist_neighbhors(h_max);
  br::smooth_hist_neighbhors(h_rate);
  br::smooth_hist_horizontal(h_avg);
  br::smooth_hist_horizontal(h_max);
  br::smooth_hist_horizontal(h_rate);

  auto* hist_dir = dir->mkdir("hists");
  hist_dir->cd();
  h_avg->Write();
  h_max->Write();
  h_rate->Write();
  h_avg_r->Write();
  h_nchan_r->Write();
  layer_accum.h_hits_evt->Write();

  draw_and_write(dir, h_avg, "c_avg", true, false);
  draw_and_write(dir, h_max, "c_max", true, false);
  draw_and_write(dir, h_rate, "c_rate", true, false);
  draw_and_write(dir, h_avg_r, "c_avg_r", false, false);
  draw_and_write(dir, layer_accum.h_hits_evt, "c_hits_evt", false, true);
}

void write_reco(TFile& output,
                const std::vector<br::occupancy::LayerAccum>& reco_layers,
                std::uint64_t n_events) {
  auto* parent = output.mkdir("reco");
  for (int layer = 0; layer <= br::occupancy::kNLayers; ++layer) {
    write_layer_directory(parent, reco_layers[layer], layer, n_events);
  }
}

void write_truth(TFile& output,
                 const std::vector<br::occupancy::TruthOccupancyGroup>& truth_groups,
                 std::uint64_t n_events) {
  for (const auto& truth_group : truth_groups) {
    auto* parent = output.mkdir(truth_group.label.c_str());
    for (int layer = 0; layer <= br::occupancy::kNLayers; ++layer) {
      write_layer_directory(parent, truth_group.layers[layer], layer, n_events);
    }
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and initialize pointers.
  const auto args = parse_args(argc, argv);
  const auto& input_dir = args.input_dir;
  const auto& output_file = args.output_file;
  const double threshold_geV = args.threshold_geV;

  // Create a decoder for cell IDs.
  const br::LFHCALDecoder decoder;
  // Find input files.
  const auto files = br::find_root_files(input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << input_dir << "\n";
    return 1;
  }

  // Ensure output directory.
  fs::path output_path = output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Initialize layer accumulators and histograms.
  std::vector<br::occupancy::LayerAccum> reco_layers;
  br::occupancy::init_reco_layers(reco_layers);
  std::vector<br::occupancy::TruthOccupancyGroup> truth_groups;
  br::occupancy::init_truth_groups(truth_groups);

  std::uint64_t n_reco_events = 0;
  std::uint64_t n_truth_events = 0;

  br::FileProgress progress(files.size(), std::cerr);

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

      if (br::occupancy::process_reco_event(frame, decoder, threshold_geV, reco_layers)) ++n_reco_events;
      if (br::occupancy::process_truth_event(frame, decoder, threshold_geV, truth_groups)) ++n_truth_events;
    }
  }

  if (n_reco_events == 0 && n_truth_events == 0) {
    std::cerr << "No events with LFHCALRecHits or LFHCALHits found\n";
    return 1;
  }

  TFile output(output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << output_file << "\n";
    return 1;
  }

  if (n_truth_events > 0) write_truth(output, truth_groups, n_truth_events);
  if (n_reco_events > 0) write_reco(output, reco_layers, n_reco_events);

  std::cout << "\n";
  output.Close();
  return 0;
}
