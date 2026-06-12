/*

./build/longitudinal_occupancy -i data/reco_bkg_apr -o plots/occupancy/longitudinal_occupancy.root

*/

#include <TCanvas.h>
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
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kRecoHitCollection = "LFHCALRecHits";
constexpr const char* kTruthHitCollection = "LFHCALHits";
constexpr double kDefaultThresholdGeV = 0.0;
constexpr double kDefaultXHalfWidthMM = 25.0;
constexpr double kEventWindowSec = 2e-6;
constexpr int kNLayers = 7;
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
  double x_half_width_mm = kDefaultXHalfWidthMM;
};

struct ChannelStats {
  double y_mm = 0.0;
  bool has_position = false;
  std::uint64_t total_hits = 0;
  int max_hits_event = 0;

  void set_y(double y) {
    y_mm = y;
    has_position = true;
  }
};

struct LayerAccum {
  std::unordered_map<br::LFHCALChannelID, ChannelStats, br::LFHCALChannelIDHash> channels;
  TH1D* h_hits_evt = nullptr;
};

struct TruthGroup {
  std::string label;
  std::vector<LayerAccum> layers;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0
            << " -i INPUT_DIR -o OUTPUT.root [-t THRESHOLD_GEV] [-x X_HALF_WIDTH_MM]\n";
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
    } else if ((arg == "-x" || arg == "--x-half-width-mm") && i + 1 < argc) {
      args.x_half_width_mm = std::stod(argv[++i]);
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

std::vector<double> make_y_edges(const LayerAccum& layer_accum, double half_width_mm) {
  std::set<double> y_coords;
  for (const auto& [channel_id, stats] : layer_accum.channels) {
    (void)channel_id;
    if (stats.has_position) y_coords.insert(stats.y_mm);
  }

  if (y_coords.empty()) return {-half_width_mm, half_width_mm};

  std::vector<double> values(y_coords.begin(), y_coords.end());
  if (values.size() == 1) return {values.front() - half_width_mm, values.front() + half_width_mm};

  std::vector<double> edges;
  edges.reserve(values.size() + 1);
  edges.push_back(values.front() - 0.5 * (values[1] - values[0]));
  for (std::size_t i = 0; i + 1 < values.size(); ++i) {
    edges.push_back(0.5 * (values[i] + values[i + 1]));
  }
  edges.push_back(values.back() + 0.5 * (values.back() - values[values.size() - 2]));
  return edges;
}

std::vector<double> make_layer_edges() {
  std::vector<double> edges(kNLayers + 1);
  for (int i = 0; i <= kNLayers; ++i) edges[i] = -0.5 + static_cast<double>(i);
  return edges;
}

void init_layers(std::vector<LayerAccum>& layers) {
  layers.assign(kNLayers, {});
  for (int layer = 0; layer < kNLayers; ++layer) {
    layers[layer].h_hits_evt = new TH1D(
        "h_hits_evt",
        (std::string("LFHCAL longitudinal slice layer ") + std::to_string(layer) + ";hits/event;Events").c_str(),
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

bool in_central_slice(const br::LFHCALCellPosition& position, double x_half_width_mm) {
  return std::abs(position.x_mm) <= x_half_width_mm;
}

bool process_reco_event(const podio::Frame& frame,
                        const br::LFHCALCellIDDecoder& decoder,
                        double threshold_geV,
                        double x_half_width_mm,
                        std::vector<LayerAccum>& layers) {
  if (!br::has_collection(frame, kRecoHitCollection)) return false;

  std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>> event_counts(kNLayers);
  std::vector<int> layer_totals(kNLayers, 0);

  const auto& hits = frame.get<edm4eic::CalorimeterHitCollection>(kRecoHitCollection);
  for (const auto& hit : hits) {
    if (hit.getEnergy() <= threshold_geV) continue;

    const int layer = hit.getLayer();
    if (layer < 0 || layer >= kNLayers) continue;

    const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
    const auto position = decoder.position(cell_id);
    if (!in_central_slice(position, x_half_width_mm)) continue;

    const auto channel_id = decoder.channel(cell_id);
    auto& stats = layers[layer].channels[channel_id];
    stats.set_y(position.y_mm);
    ++event_counts[layer][channel_id];
    ++layer_totals[layer];
  }

  for (int layer = 0; layer < kNLayers; ++layer) {
    layers[layer].h_hits_evt->Fill(layer_totals[layer]);
    for (const auto& [channel_id, count] : event_counts[layer]) {
      auto& stats = layers[layer].channels[channel_id];
      stats.total_hits += static_cast<std::uint64_t>(count);
      stats.max_hits_event = std::max(stats.max_hits_event, count);
    }
  }

  return true;
}

bool process_truth_event(const podio::Frame& frame,
                         const br::LFHCALCellIDDecoder& decoder,
                         double threshold_geV,
                         double x_half_width_mm,
                         std::vector<TruthGroup>& groups) {
  if (!br::has_collection(frame, kTruthHitCollection)) return false;

  std::vector<std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>>> event_counts(
      groups.size(), std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>>(kNLayers));
  std::vector<std::vector<int>> layer_totals(groups.size(), std::vector<int>(kNLayers, 0));

  const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kTruthHitCollection);
  for (const auto& hit : hits) {
    if (hit.getEnergy() <= threshold_geV) continue;

    const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
    const auto channel_id = decoder.channel(cell_id);
    const int layer = channel_id.rlayerz;
    if (layer < 0 || layer >= kNLayers) continue;

    const auto position = decoder.position(cell_id);
    if (!in_central_slice(position, x_half_width_mm)) continue;

    std::array<double, kAllTruthIndex> energy_by_origin{};
    for (const auto& contribution : hit.getContributions()) {
      const double contribution_energy = contribution.getEnergy();
      if (contribution_energy <= 0.0) continue;
      const int origin_index = br::origin_index(contribution.getParticle().getGeneratorStatus());
      if (origin_index < 0 || origin_index >= kAllTruthIndex) continue;
      energy_by_origin[origin_index] += contribution_energy;
    }

    for (int origin_index = 0; origin_index < kAllTruthIndex; ++origin_index) {
      if (energy_by_origin[origin_index] <= 0.0) continue;
      auto& stats = groups[origin_index].layers[layer].channels[channel_id];
      stats.set_y(position.y_mm);
      ++event_counts[origin_index][layer][channel_id];
      ++layer_totals[origin_index][layer];
    }

    auto& all_stats = groups[kAllTruthIndex].layers[layer].channels[channel_id];
    all_stats.set_y(position.y_mm);
    ++event_counts[kAllTruthIndex][layer][channel_id];
    ++layer_totals[kAllTruthIndex][layer];
  }

  for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
    for (int layer = 0; layer < kNLayers; ++layer) {
      groups[group_index].layers[layer].h_hits_evt->Fill(layer_totals[group_index][layer]);
      for (const auto& [channel_id, count] : event_counts[group_index][layer]) {
        auto& stats = groups[group_index].layers[layer].channels[channel_id];
        stats.total_hits += static_cast<std::uint64_t>(count);
        stats.max_hits_event = std::max(stats.max_hits_event, count);
      }
    }
  }

  return true;
}

void write_group(TDirectory* parent,
                 const std::vector<LayerAccum>& layers,
                 const std::string& title_prefix,
                 std::uint64_t n_events,
                 double x_half_width_mm,
                 const std::vector<double>* reference_y_edges = nullptr) {
  std::set<double> all_y;
  double rmax = 0.0;
  for (const auto& layer_accum : layers) {
    for (const auto& [channel_id, stats] : layer_accum.channels) {
      (void)channel_id;
      if (stats.has_position) all_y.insert(stats.y_mm);
      rmax = std::max(rmax, std::abs(stats.y_mm));
    }
  }

  LayerAccum merged;
  for (double y : all_y) {
    br::LFHCALChannelID dummy{};
    dummy.moduleIDx = static_cast<int>(y);
    auto& stats = merged.channels[dummy];
    stats.set_y(y);
  }

  const auto y_edges = reference_y_edges ? *reference_y_edges : make_y_edges(merged, x_half_width_mm);
  const auto layer_edges = make_layer_edges();

  auto* h_avg = new TH2D(
      "h_avg",
      (title_prefix + ";readout layer;y [mm];avg hits/event/channel").c_str(),
      kNLayers, layer_edges.data(), static_cast<int>(y_edges.size()) - 1, y_edges.data());
  auto* h_max = new TH2D(
      "h_max",
      (title_prefix + ";readout layer;y [mm];max hits/event/channel").c_str(),
      kNLayers, layer_edges.data(), static_cast<int>(y_edges.size()) - 1, y_edges.data());
  auto* h_rate = new TH2D(
      "h_rate",
      (title_prefix + ";readout layer;y [mm];rate [Hz/channel]").c_str(),
      kNLayers, layer_edges.data(), static_cast<int>(y_edges.size()) - 1, y_edges.data());
  auto* h_avg_y = new TH1D(
      "h_avg_y",
      (title_prefix + ";y [mm];avg hits/event").c_str(),
      static_cast<int>(y_edges.size()) - 1, y_edges.data());
  auto* h_nchan_y = new TH1D("h_nchan_y", "", static_cast<int>(y_edges.size()) - 1, y_edges.data());

  auto* per_layer_dir = parent->mkdir("layers");
  for (int layer = 0; layer < kNLayers; ++layer) {
    const auto& layer_accum = layers[layer];
    auto* layer_dir = per_layer_dir->mkdir(("layer" + std::to_string(layer)).c_str());
    auto* hist_dir = layer_dir->mkdir("hists");
    hist_dir->cd();
    layer_accum.h_hits_evt->Write();
  }

  for (int layer = 0; layer < kNLayers; ++layer) {
    for (const auto& [channel_id, stats] : layers[layer].channels) {
      (void)channel_id;
      const double avg = static_cast<double>(stats.total_hits) / static_cast<double>(n_events);
      const double rate = static_cast<double>(stats.total_hits) / (static_cast<double>(n_events) * kEventWindowSec);
      h_avg->Fill(static_cast<double>(layer), stats.y_mm, avg);
      h_max->Fill(static_cast<double>(layer), stats.y_mm, stats.max_hits_event);
      h_rate->Fill(static_cast<double>(layer), stats.y_mm, rate);
      h_avg_y->Fill(stats.y_mm, avg);
      h_nchan_y->Fill(stats.y_mm, 1.0);
    }
  }

  for (int bin = 1; bin <= h_avg_y->GetNbinsX(); ++bin) {
    const double nch = h_nchan_y->GetBinContent(bin);
    if (nch > 0.0) h_avg_y->SetBinContent(bin, h_avg_y->GetBinContent(bin) / nch);
  }

  auto* hist_dir = parent->mkdir("hists");
  hist_dir->cd();
  h_avg->Write();
  h_max->Write();
  h_rate->Write();
  h_avg_y->Write();
  h_nchan_y->Write();

  draw_and_write(parent, h_avg, "c_avg", true, false);
  draw_and_write(parent, h_max, "c_max", true, false);
  draw_and_write(parent, h_rate, "c_rate", true, false);
  draw_and_write(parent, h_avg_y, "c_avg_y", false, false);
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  const auto args = parse_args(argc, argv);
  const auto& input_dir = args.input_dir;
  const auto& output_file = args.output_file;
  const double threshold_geV = args.threshold_geV;
  const double x_half_width_mm = args.x_half_width_mm;

  const br::LFHCALCellIDDecoder decoder;
  const auto files = br::find_root_files(input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << input_dir << "\n";
    return 1;
  }

  fs::path output_path = output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  std::vector<LayerAccum> reco_layers;
  init_layers(reco_layers);
  std::vector<TruthGroup> truth_groups;
  init_truth_groups(truth_groups);

  std::uint64_t n_reco_events = 0;
  std::uint64_t n_truth_events = 0;

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
      if (process_reco_event(frame, decoder, threshold_geV, x_half_width_mm, reco_layers)) ++n_reco_events;
      if (process_truth_event(frame, decoder, threshold_geV, x_half_width_mm, truth_groups)) ++n_truth_events;
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

  if (n_reco_events > 0) {
    auto* reco_dir = output.mkdir("reco");
    write_group(reco_dir, reco_layers, "Reco longitudinal occupancy", n_reco_events, x_half_width_mm);
  }

  if (n_truth_events > 0) {
    std::vector<double> dis_reference_y_edges;
    bool have_dis_reference = false;
    for (const auto& group : truth_groups) {
      if (group.label == "DIS") {
        std::set<double> dis_y;
        LayerAccum merged;
        for (const auto& layer_accum : group.layers) {
          for (const auto& [channel_id, stats] : layer_accum.channels) {
            (void)channel_id;
            if (stats.has_position) dis_y.insert(stats.y_mm);
          }
        }
        for (double y : dis_y) {
          br::LFHCALChannelID dummy{};
          dummy.moduleIDx = static_cast<int>(y);
          auto& stats = merged.channels[dummy];
          stats.set_y(y);
        }
        dis_reference_y_edges = make_y_edges(merged, x_half_width_mm);
        have_dis_reference = true;
        break;
      }
    }

    for (const auto& group : truth_groups) {
      auto* dir = output.mkdir(group.label.c_str());
      write_group(dir,
                  group.layers,
                  group.label + " longitudinal occupancy",
                  n_truth_events,
                  x_half_width_mm,
                  have_dis_reference ? &dis_reference_y_edges : nullptr);
    }
  }

  std::cout << "\n";
  output.Close();
  return 0;
}
