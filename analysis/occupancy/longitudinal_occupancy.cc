/*

./build/longitudinal_occupancy -i data/bkg_apr -o plots/occupancy/longitudinal_occupancy.root

*/

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include "truth_longitudinal.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
namespace lo = br::occupancy::longitudinal;

namespace {

// Constant(s)
constexpr double MIP_1 = 3.5e-3;
constexpr double MIP_2 = 7.0e-3;
constexpr double kCoefficient = 0.5;
const lo::ThresholdsByLayer kThresholdsGeV = {
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

struct YBinStats {
  double sum_avg = 0.0;
  double y_mm = 0.0;
  int n_channels = 0;
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

void draw_and_write(TDirectory* dir, TH1* hist, const char* canvas_name, bool logz = false, bool logy = false) {
  dir->cd();
  TCanvas canvas(canvas_name, hist->GetTitle(), 1000, 800);
  if (logz) canvas.SetLogz();
  if (logy) canvas.SetLogy();
  hist->SetStats(false);
  hist->Draw(hist->InheritsFrom(TH2::Class()) ? "colz" : "hist");
  canvas.Write();
}

std::vector<double> make_axis_edges(const std::set<double>& coords) {
  // Build ROOT bin edges from channel center coordinates.
  if (coords.empty()) return {-lo::kDisplayPaddingMM, lo::kDisplayPaddingMM};

  std::vector<double> values(coords.begin(), coords.end());

  // With a single row, make one padded bin around the channel center.
  if (values.size() == 1) {
    return {values.front() - lo::kDisplayPaddingMM, values.front() + lo::kDisplayPaddingMM};
  }

  // Use the smallest observed channel spacing as the detector row pitch.
  // Larger jumps are geometry gaps, such as the beamline opening, so they
  // should become empty bins in the final histogram.
  double min_step = values[1] - values[0];
  for (std::size_t i = 1; i < values.size(); ++i) {
    min_step = std::min(min_step, values[i] - values[i - 1]);
  }

  // Put each observed channel center in the middle of a pitch-sized bin.
  const double half_step = 0.5 * min_step;
  const double first_edge = values.front() - half_step;
  const double last_edge = values.back() + half_step;
  const int nbins = static_cast<int>(std::lround((last_edge - first_edge) / min_step));

  // Fill every pitch step between the extremes. Missing channel rows are
  // represented by empty bins, which leaves visual gaps in the y direction.
  std::vector<double> edges;
  edges.reserve(nbins + 1);
  for (int i = 0; i <= nbins; ++i) {
    edges.push_back(first_edge + static_cast<double>(i) * min_step);
  }
  return edges;
}

std::vector<double> collect_y_edges(const std::vector<lo::LayerAccum>& layers) {
  // Collect every y row seen in any readout layer for this origin group.
  // Sharing one y axis across layers makes the y-z view line up vertically.
  std::set<double> y_coords;
  for (const auto& layer : layers) {
    for (const auto& [channel_id, stats] : layer.channels) {
      (void)channel_id;
      y_coords.insert(stats.y_mm);
    }
  }
  return make_axis_edges(y_coords);
}

void write_group_directory(TDirectory* parent,
                           const std::vector<lo::LayerAccum>& layers,
                           std::uint64_t n_events,
                           const std::optional<std::vector<double>>& reference_y_edges = std::nullopt) {
  // Build common y binning across layers for this group.
  std::vector<double> y_edges = reference_y_edges.has_value() ? *reference_y_edges : collect_y_edges(layers);

  // These are the longitudinal equivalents of the transverse occupancy products:
  // layer on x, y on y, and channel-averaged content in each (layer, y) bin.
  auto* h_avg = new TH2D(
      "h_avg",
      "Occupancy in central slice;readout layer;y [mm];avg hits/event/channel",
      lo::kNLayers, -0.5, lo::kNLayers - 0.5,
      static_cast<int>(y_edges.size()) - 1, y_edges.data());

  auto* h_rate = new TH2D(
      "h_rate",
      "Rates in central slice;readout layer;y [mm];rate [Hz/channel]",
      lo::kNLayers, -0.5, lo::kNLayers - 0.5,
      static_cast<int>(y_edges.size()) - 1, y_edges.data());

  auto* h_avg_y = new TH1D(
      "h_avg_y",
      "Occupancy in central slice;y [mm];avg hits/event/channel",
      static_cast<int>(y_edges.size()) - 1, y_edges.data());

  // One map per readout layer; each map groups channels by y row.
  std::vector<std::unordered_map<double, YBinStats>> yz_bins(lo::kNLayers);
  std::unordered_map<double, YBinStats> y_bins;

  // First collect all channels that land in the same (layer, y) bin,
  // so channels within the x slice get averaged together.
  for (int layer = 0; layer < lo::kNLayers; ++layer) {
    for (const auto& [channel_id, stats] : layers[layer].channels) {
      (void)channel_id;
      const double avg = static_cast<double>(stats.total_hits) / static_cast<double>(n_events);

      // yz bin is defined by the readout layer and the y position of the channel.
      auto& yz = yz_bins[layer][stats.y_mm];
      yz.y_mm = stats.y_mm;
      yz.sum_avg += avg;
      ++yz.n_channels;

      auto& y = y_bins[stats.y_mm];
      y.y_mm = stats.y_mm;
      y.sum_avg += avg;
      ++y.n_channels;
    }
  }

  // Then write one averaged value per (layer, y) bin into the 2D histograms.
  for (int layer = 0; layer < lo::kNLayers; ++layer) {
    for (const auto& [y_mm, stats] : yz_bins[layer]) {
      (void)y_mm;
      if (stats.n_channels == 0) continue;
      const double avg = stats.sum_avg / static_cast<double>(stats.n_channels);
      h_avg->Fill(layer, stats.y_mm, avg);
      h_rate->Fill(layer, stats.y_mm, avg / lo::kEventWindowSec);
    }
  }

  // Also collapse over z to make a y projection of the average occupancy.
  for (const auto& [y_mm, stats] : y_bins) {
    (void)y_mm;
    if (stats.n_channels == 0) continue;
    h_avg_y->Fill(stats.y_mm, stats.sum_avg / static_cast<double>(stats.n_channels));
  }

  auto* hist_dir = parent->mkdir("hists");
  hist_dir->cd();
  h_avg->Write();
  h_rate->Write();
  h_avg_y->Write();
  for (int layer = 0; layer < lo::kNLayers; ++layer) {
    layers[layer].h_hits_evt->Write();
  }

  draw_and_write(parent, h_avg, "c_avg", true, false);
  draw_and_write(parent, h_rate, "c_rate", true, false);
  draw_and_write(parent, h_avg_y, "c_avg_y", false, false);
  for (int layer = 0; layer < lo::kNLayers; ++layer) {
    draw_and_write(parent, layers[layer].h_hits_evt, ("c_hits_evt_layer" + std::to_string(layer)).c_str(), false, true);
  }
}

void write_truth(TFile& output,
                 const std::vector<lo::TruthGroup>& truth_groups,
                 std::uint64_t n_events) {
  std::optional<std::vector<double>> dis_y_edges = std::nullopt;
  const auto dis_it = std::find_if(
      truth_groups.begin(),
      truth_groups.end(),
      [](const lo::TruthGroup& group) { return group.label == "DIS"; });
  if (dis_it != truth_groups.end()) dis_y_edges = collect_y_edges(dis_it->layers);

  for (const auto& group : truth_groups) {
    auto* dir = output.mkdir(group.label.c_str());
    write_group_directory(dir, group.layers, n_events, dis_y_edges);
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and gather input files.
  const auto args = parse_args(argc, argv);
  const auto files = br::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Create a decoder for cell IDs.
  const br::LFHCALCellIDDecoder decoder;

  // NOTE: reco unavailable for now due to some timing issues with background.
  std::vector<lo::TruthGroup> truth_groups;
  lo::init_truth_groups(truth_groups);

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

      // NOTE: reco unavailable for now due to some timing issues with background.
      if (lo::process_truth_event(frame, decoder, kThresholdsGeV, truth_groups)) ++n_truth_events;
    }
  }

  if (n_truth_events == 0) {
    std::cerr << "No events with LFHCALHits found\n";
    return 1;
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  // Write output plots.
  if (n_truth_events > 0) write_truth(output, truth_groups, n_truth_events);

  std::cout << "\n";
  output.Close();
  return 0;
}
