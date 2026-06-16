/*

./build/longitudinal_occupancy -i data/reco_bkg_apr -o plots/occupancy/longitudinal_occupancy.root

*/

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4eic/CalorimeterHitCollection.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "classify_hit.h"
#include "decode_cell_id.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kDefaultThresholdGeV = 0.001;
constexpr int kNLayers = 7;
constexpr double kEventWindowSec = 2e-6;
constexpr double kDisplayPaddingMM = 25.0;
constexpr double kXMinMM = -400.0;
constexpr double kXMaxMM = 300.0;
constexpr const char* kRecoHitCollection = "LFHCALRecHits";
constexpr const char* kTruthHitCollection = "LFHCALHits";

constexpr const char* kTruthLabels[] = {
    "DIS",
    "synrad",
    "eBrem",
    "eTouschek",
    "eCoulomb",
    "pBeamGas",
    "other",
    "all",
};

constexpr int kNTruthGroups = sizeof(kTruthLabels) / sizeof(kTruthLabels[0]);
constexpr int kAllTruthIndex = kNTruthGroups - 1;

struct Args {
  std::string input_dir;
  std::string output_file;
  double threshold_geV = kDefaultThresholdGeV;
};

struct ChannelStats {
  double y_mm = 0.0;
  std::uint64_t total_hits = 0;
  int max_hits_event = 0;
};

struct LayerAccum {
  std::unordered_map<br::LFHCALChannelID, ChannelStats, br::LFHCALChannelIDHash> channels;
  TH1D* h_hits_evt = nullptr;
};

struct TruthGroup {
  std::string label;
  std::vector<LayerAccum> layers;
};

struct YBinStats {
  double sum_avg = 0.0;
  double sum_max = 0.0;
  double y_mm = 0.0;
  int n_channels = 0;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT_DIR -o OUTPUT.root [-t THRESHOLD_GEV]\n";
}

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
  if (coords.empty()) return {-kDisplayPaddingMM, kDisplayPaddingMM};

  std::vector<double> values(coords.begin(), coords.end());

  // With a single row, make one padded bin around the channel center.
  if (values.size() == 1) {
    return {values.front() - kDisplayPaddingMM, values.front() + kDisplayPaddingMM};
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
  // represented by empty bins, which should leaves visual gaps in the y direction.
  std::vector<double> edges;
  edges.reserve(nbins + 1);
  for (int i = 0; i <= nbins; ++i) {
    edges.push_back(first_edge + static_cast<double>(i) * min_step);
  }
  return edges;
}

bool in_x_slice(const br::LFHCALCellPosition& position) {
  return position.x_mm >= kXMinMM && position.x_mm < kXMaxMM;
}

void init_layers(std::vector<LayerAccum>& layers) {
  layers.assign(kNLayers, {});
  for (int layer = 0; layer < kNLayers; ++layer) {
    layers[layer].h_hits_evt = new TH1D(
        ("h_hits_evt_layer" + std::to_string(layer)).c_str(),
        ("LFHCAL layer " + std::to_string(layer) + ";hits/event;Events").c_str(),
        200, 0, 200);
  }
}

void init_truth_groups(std::vector<TruthGroup>& groups) {
  groups.clear();
  groups.reserve(kNTruthGroups);
  for (const char* label : kTruthLabels) {
    TruthGroup group;
    group.label = label;
    init_layers(group.layers);
    groups.push_back(std::move(group));
  }
}

bool process_reco_event(const podio::Frame& frame,
                        const br::LFHCALCellIDDecoder& decoder,
                        double threshold_geV,
                        std::vector<LayerAccum>& layers) {
  if (!br::has_collection(frame, kRecoHitCollection)) return false;

  // Per-event channel counters.
  std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>> event_counts(kNLayers);
  std::vector<int> layer_totals(kNLayers, 0);

  // Grab and loop over hits.
  const auto& hits = frame.get<edm4eic::CalorimeterHitCollection>(kRecoHitCollection);
  for (const auto& hit : hits) {
    // Apply threshold.
    if (hit.getEnergy() <= threshold_geV) continue;

    // Get layer.
    const int layer = hit.getLayer();
    if (layer < 0 || layer >= kNLayers) continue;

    // Decode cell ID and apply the x slice.
    const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
    const auto position = decoder.position(cell_id);
    if (!in_x_slice(position)) continue;

    // All channels with the same y value are stacked on top of each other.
    const auto channel_id = decoder.channel(cell_id);
    auto& stats = layers[layer].channels[channel_id];
    stats.y_mm = position.y_mm;

    // Count one more hit in this channel for the current event.
    ++event_counts[layer][channel_id];
    ++layer_totals[layer];
  }

  // End of this event:
  // move the temporary per-event counts into the long-lived channel statistics.
  for (int layer = 0; layer < kNLayers; ++layer) {
    layers[layer].h_hits_evt->Fill(layer_totals[layer]);
    for (const auto& [channel_id, count] : event_counts[layer]) {
      auto& stats = layers[layer].channels[channel_id];

      // Add this event's contribution to the all-events total.
      stats.total_hits += static_cast<std::uint64_t>(count);

      // Keep the largest number of hits this channel ever saw in any single event.
      stats.max_hits_event = std::max(stats.max_hits_event, count);
    }
  }

  return true;
}

bool process_truth_event(const podio::Frame& frame,
                         const br::LFHCALCellIDDecoder& decoder,
                         double threshold_geV,
                         std::vector<TruthGroup>& groups) {
  // Try to grab LFHCALHits.
  if (!br::has_collection(frame, kTruthHitCollection)) return false;

  // Per-event channel counters for every truth group and layer.
  std::vector<std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>>> event_counts(
      groups.size(), std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>>(kNLayers));
  std::vector<std::vector<int>> layer_totals(groups.size(), std::vector<int>(kNLayers, 0));

  // Grab and loop over hits.
  const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kTruthHitCollection);
  for (const auto& hit : hits) {
    // Apply threshold.
    if (hit.getEnergy() <= threshold_geV) continue;

    // Decode cell ID and apply the x slice.
    const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
    const auto position = decoder.position(cell_id);
    if (!in_x_slice(position)) continue;

    // Get readout layer and channel ID for this hit.
    const auto channel_id = decoder.channel(cell_id);
    const int layer = channel_id.rlayerz;
    if (layer < 0 || layer >= kNLayers) continue;

    std::array<double, kAllTruthIndex> energy_by_origin{};

    // Per hit, sum all contributions that fall into the same broad origin family.
    for (const auto& contribution : hit.getContributions()) {
      const double energy = contribution.getEnergy();
      if (energy <= 0.0) continue;

      const int origin_index = br::origin_index(contribution.getParticle().getGeneratorStatus());
      if (origin_index < 0 || origin_index >= kAllTruthIndex) continue;
      energy_by_origin[origin_index] += energy;
    }

    // Loop through all seen origins for this hit.
    for (int origin_index = 0; origin_index < kAllTruthIndex; ++origin_index) {
      // Skip if this origin contributes no energy.
      if (energy_by_origin[origin_index] <= 0.0) continue;

      // Otherwise fill the corresponding longitudinal bucket for this origin.
      auto& stats = groups[origin_index].layers[layer].channels[channel_id];
      stats.y_mm = position.y_mm;
      ++event_counts[origin_index][layer][channel_id];
      ++layer_totals[origin_index][layer];
    }

    // Also update the inclusive truth bucket once per hit.
    auto& all_stats = groups[kAllTruthIndex].layers[layer].channels[channel_id];
    all_stats.y_mm = position.y_mm;
    ++event_counts[kAllTruthIndex][layer][channel_id];
    ++layer_totals[kAllTruthIndex][layer];
  }

  // For each generatorStatus family...
  for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
    // For each readout layer...
    for (int layer = 0; layer < kNLayers; ++layer) {
      groups[group_index].layers[layer].h_hits_evt->Fill(layer_totals[group_index][layer]);

      // For each channel, accumulate stats.
      for (const auto& [channel_id, count] : event_counts[group_index][layer]) {
        auto& stats = groups[group_index].layers[layer].channels[channel_id];

        // Add this event's contribution to the all-events total.
        stats.total_hits += static_cast<std::uint64_t>(count);

        // Keep track of the single-event maximum for this channel.
        stats.max_hits_event = std::max(stats.max_hits_event, count);
      }
    }
  }

  return true;
}

std::vector<double> collect_y_edges(const std::vector<LayerAccum>& layers) {
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
                           const std::vector<LayerAccum>& layers,
                           std::uint64_t n_events,
                           const std::optional<std::vector<double>>& reference_y_edges = std::nullopt) {
  // Build common y binning across layers for this group.
  std::vector<double> y_edges = reference_y_edges.has_value() ? *reference_y_edges : collect_y_edges(layers);

  // These are the longitudinal equivalents of the transverse occupancy products:
  // layer on x, y on y, and channel-averaged content in each (layer, y) bin.
  auto* h_avg = new TH2D(
      "h_avg",
      "Longitudinal occupancy;readout layer;y [mm];avg hits/event/channel",
      kNLayers, -0.5, kNLayers - 0.5,
      static_cast<int>(y_edges.size()) - 1, y_edges.data());

  auto* h_max = new TH2D(
      "h_max",
      "Longitudinal occupancy;readout layer;y [mm];avg max hits/event/channel",
      kNLayers, -0.5, kNLayers - 0.5,
      static_cast<int>(y_edges.size()) - 1, y_edges.data());

  auto* h_rate = new TH2D(
      "h_rate",
      "Longitudinal occupancy;readout layer;y [mm];rate [Hz/channel]",
      kNLayers, -0.5, kNLayers - 0.5,
      static_cast<int>(y_edges.size()) - 1, y_edges.data());

  auto* h_avg_y = new TH1D(
      "h_avg_y",
      "Longitudinal occupancy;y [mm];avg hits/event/channel",
      static_cast<int>(y_edges.size()) - 1, y_edges.data());

  // One map per readout layer; each map groups channels by y row.
  std::vector<std::unordered_map<double, YBinStats>> yz_bins(kNLayers);
  std::unordered_map<double, YBinStats> y_bins;

  // First collect all channels that land in the same (layer, y) bin,
  // so channels within the x slice get averaged together.
  for (int layer = 0; layer < kNLayers; ++layer) {
    for (const auto& [channel_id, stats] : layers[layer].channels) {
      (void)channel_id;
      const double avg = static_cast<double>(stats.total_hits) / static_cast<double>(n_events);

      // yz bin is defined by the readout layer and the y position of the channel.
      auto& yz = yz_bins[layer][stats.y_mm];
      yz.y_mm = stats.y_mm;
      yz.sum_avg += avg;
      yz.sum_max += static_cast<double>(stats.max_hits_event);
      ++yz.n_channels;

      auto& y = y_bins[stats.y_mm];
      y.y_mm = stats.y_mm;
      y.sum_avg += avg;
      ++y.n_channels;
    }
  }

  // Then write one averaged value per (layer, y) bin into the 2D histograms.
  for (int layer = 0; layer < kNLayers; ++layer) {
    for (const auto& [y_mm, stats] : yz_bins[layer]) {
      (void)y_mm;
      if (stats.n_channels == 0) continue;
      const double avg = stats.sum_avg / static_cast<double>(stats.n_channels);
      const double avg_max = stats.sum_max / static_cast<double>(stats.n_channels);
      h_avg->Fill(layer, stats.y_mm, avg);
      h_max->Fill(layer, stats.y_mm, avg_max);
      h_rate->Fill(layer, stats.y_mm, avg / kEventWindowSec);
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
  h_max->Write();
  h_rate->Write();
  h_avg_y->Write();
  for (int layer = 0; layer < kNLayers; ++layer) {
    layers[layer].h_hits_evt->Write();
    draw_and_write(parent, layers[layer].h_hits_evt, ("c_hits_evt_layer" + std::to_string(layer)).c_str(), false, true);
  }

  draw_and_write(parent, h_avg, "c_avg", true, false);
  draw_and_write(parent, h_max, "c_max", true, false);
  draw_and_write(parent, h_rate, "c_rate", true, false);
  draw_and_write(parent, h_avg_y, "c_avg_y", false, false);
}

void write_reco(TFile& output,
                const std::vector<LayerAccum>& reco_layers,
                std::uint64_t n_events) {
  auto* dir = output.mkdir("reco");
  write_group_directory(dir, reco_layers, n_events);
}

void write_truth(TFile& output,
                 const std::vector<TruthGroup>& truth_groups,
                 std::uint64_t n_events) {
  std::optional<std::vector<double>> dis_y_edges = std::nullopt;
  const auto dis_it = std::find_if(
      truth_groups.begin(),
      truth_groups.end(),
      [](const TruthGroup& group) { return group.label == "DIS"; });
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

  // Initialize accumulators for reco and truth.
  std::vector<LayerAccum> reco_layers;
  init_layers(reco_layers);
  std::vector<TruthGroup> truth_groups;
  init_truth_groups(truth_groups);

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

      // Process reco and truth occupancy in the longitudinal view.
      if (process_reco_event(frame, decoder, args.threshold_geV, reco_layers)) ++n_reco_events;
      if (process_truth_event(frame, decoder, args.threshold_geV, truth_groups)) ++n_truth_events;
    }
  }

  if (n_reco_events == 0 && n_truth_events == 0) {
    std::cerr << "No events with LFHCALRecHits or LFHCALHits found\n";
    return 1;
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  // Write output plots.
  if (n_truth_events > 0) write_truth(output, truth_groups, n_truth_events);
  if (n_reco_events > 0) write_reco(output, reco_layers, n_reco_events);

  std::cout << "\n";
  output.Close();
  return 0;
}
