/*

./build/occupancy -i data/reco -o plots/occupancy/occupancy.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4eic/CalorimeterHitCollection.h>

#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

// Constants
constexpr const char* kHitCollection = "LFHCALRecHits"; // Use reco collection
constexpr double kDefaultThresholdGeV = 5e-4;           // Default threshold is 0.5 MeV.
constexpr double kEventWindowSec = 2e-6;                // 2 microseconds.
constexpr int kNLayers = 7;                             // 7 longitudinal readout layers.
constexpr int kAllLayersIndex = kNLayers;               // Index for merged layer stats.
constexpr double kCellSizeMM = 50.0;                    // Cell size is 5cm by 5 cm.
constexpr double kXYExtentMM = 2700.0;                  // +- 2700 mm for x and y ranges.

// Stats per readout layer.
struct ChannelStats {
  double x = 0.0;
  double y = 0.0;
  double r = 0.0;
  std::uint64_t total_hits = 0;
  int max_hits_event = 0;
};

// Stores stats for all 7 layers plus a summed layer consisting of all 7 layers.
struct LayerAccum {
  // For one layer bucket, map each unique transverse cellID to accumulated channel stats.
  // Since the enclosing vector index already selects the layer, the full readout channel is
  // effectively (layer index, cellID).
  std::unordered_map<std::uint64_t, ChannelStats> channels;
  TH1D* h_hits_evt = nullptr;
};

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

// Draw a (custom) histogram on a canvas and write it to the given directory.
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

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and initialize pointers.
  const auto args = parse_args(argc, argv);
  const auto& input_dir = args.input_dir;
  const auto& output_file = args.output_file;
  const double threshold_geV = args.threshold_geV;

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
  std::vector<LayerAccum> layers(kNLayers + 1); // +1 to account for the summed layer stats.
  for (int layer = 0; layer <= kNLayers; ++layer) {
    layers[layer].h_hits_evt = new TH1D(
        "h_hits_evt",
        (layer == kAllLayersIndex ? std::string("LFHCAL summed layers;hits/event;Events")
          : std::string("LFHCAL layer ") + std::to_string(layer) + ";hits/event;Events").c_str(),
        200, 0, 200);
  }

  std::uint64_t n_events = 0;

  const auto total_files = files.size();
  std::size_t i = 0;

  // Loop over all files in the input directory.
  for (const auto& path : files) {
    ++i;
    std::cerr << "\r" << i << "/" << total_files << " files read." << std::flush;

    podio::ROOTReader reader;
    reader.openFile(path.string());

    const std::size_t total_events = reader.getEntries("events");

    // Per file, loop over all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {

      // Read event.
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      // Grab LFHCALRecHits when possible.
      podio::Frame frame(std::move(data));
      if (!br::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      // Initialize per-event counts.
      std::vector<std::unordered_map<std::uint64_t, int>> event_counts(kNLayers + 1);
      std::vector<int> layer_totals(kNLayers + 1, 0);

      // Grab and loop over hits.
      const auto& hits = frame.get<edm4eic::CalorimeterHitCollection>(kHitCollection);
      for (const auto& hit : hits) {

        // Apply hit threshold.
        if (hit.getEnergy() <= threshold_geV) continue;

        // Get layer and cellID.
        const int layer = hit.getLayer();
        if (layer < 0 || layer >= kNLayers) continue;

        // Copy and accumulate hit info to the corresponding layer channel.
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        const auto pos = hit.getPosition();

        // layers[layer].channels stores all unique cellIDs seen for that layer.
        auto& stats = layers[layer].channels[cell_id]; // So one "channel" = all cells with the same (x, y, layer).
        stats.x = pos.x;
        stats.y = pos.y;
        stats.r = std::hypot(pos.x, pos.y);

        ++event_counts[layer][cell_id];
        ++layer_totals[layer];

        // Also accumulate stats for the summed layer.
        // A synthetic key from (layer, cell_id) is built before storing it in the summed bucket.
        const std::uint64_t merged_cell_id = static_cast<std::uint64_t>(layer) * (1ULL << 56) | cell_id;
        auto& merged_stats = layers[kAllLayersIndex].channels[merged_cell_id];
        merged_stats.x = pos.x;
        merged_stats.y = pos.y;
        merged_stats.r = std::hypot(pos.x, pos.y);
        ++event_counts[kAllLayersIndex][merged_cell_id];
        ++layer_totals[kAllLayersIndex];
      }

      // Fill hists with accumulated stats.
      for (int layer = 0; layer <= kNLayers; ++layer) {
        layers[layer].h_hits_evt->Fill(layer_totals[layer]);
        for (const auto& [cell_id, count] : event_counts[layer]) {
          auto& stats = layers[layer].channels[cell_id];
          stats.total_hits += static_cast<std::uint64_t>(count);
          stats.max_hits_event = std::max(stats.max_hits_event, count);
        }
      }
    }
  }

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  TFile output(output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << output_file << "\n";
    return 1;
  }

  // Loop over layers and create histograms + canvases.
  for (int layer = 0; layer <= kNLayers; ++layer) {

    // Create directory.
    const std::string dir_name = layer == kAllLayersIndex ? "sum_layers" : std::string("layer") + std::to_string(layer);
    auto* dir = output.mkdir(dir_name.c_str());
    dir->cd();

    // If no channels has hits, just write the hits/event histogram and skip the rest.
    if (layers[layer].channels.empty()) {
      auto* hist_dir = dir->mkdir("hists");
      hist_dir->cd();
      layers[layer].h_hits_evt->Write();
      continue;
    }

    // Handle axis ranges.
    double rmax = 0.0;
    for (const auto& [cell_id, stats] : layers[layer].channels) {
      (void)cell_id;
      rmax = std::max(rmax, stats.r);
    }
    const double xmin = -kXYExtentMM;
    const double xmax = kXYExtentMM;
    const double ymin = -kXYExtentMM;
    const double ymax = kXYExtentMM;
    const int xbins = static_cast<int>(std::ceil((xmax - xmin) / kCellSizeMM));
    const int ybins = static_cast<int>(std::ceil((ymax - ymin) / kCellSizeMM));

    // -------------------------- initialize TH2Ds --------------------------
    // Average hits/event/channel.
    auto* h_avg = new TH2D(
        "h_avg",
        (layer == kAllLayersIndex ? std::string("Summed layers;x [mm];y [mm];avg hits/event/channel")
          : std::string("Layer ") + std::to_string(layer) + ";x [mm];y [mm];avg hits/event/channel").c_str(),
        xbins, xmin, xmax,
        ybins, ymin, ymax);

    // Max hits/event/channel.
    auto* h_max = new TH2D(
        "h_max",
        (layer == kAllLayersIndex ? std::string("Summed layers;x [mm];y [mm];max hits/event/channel")
          : std::string("Layer ") + std::to_string(layer) + ";x [mm];y [mm];max hits/event/channel").c_str(),
        xbins, xmin, xmax,
        ybins, ymin, ymax);

    // Rate in Hz/channel = total hits / (total events * event time window).
    auto* h_rate = new TH2D(
        "h_rate",
        (layer == kAllLayersIndex ? std::string("Summed layers;x [mm];y [mm];rate [Hz/channel]")
          : std::string("Layer ") + std::to_string(layer) + ";x [mm];y [mm];rate [Hz/channel]").c_str(),
        xbins, xmin, xmax,
        ybins, ymin, ymax);

    // Average hits/events/channel vs radius.
    auto* h_avg_r = new TH1D(
        "h_avg_r",
        (layer == kAllLayersIndex ? std::string("Summed layers;R [mm];avg hits/event")
          : std::string("Layer ") + std::to_string(layer) + ";R [mm];avg hits/event").c_str(),
        100, 0.0, std::max(1.0, 1.05 * rmax));

    // Number of channels vs radius (for normalization).
    auto* h_nchan_r = new TH1D("h_nchan_r", "", 100, 0.0, std::max(1.0, 1.05 * rmax));

    // Per channel, compute stats and fill corresponding histograms.
    // Loop over all unique channels collected for this layer bucket.
    // layer0..layer6 --> all cells that share (x, y, layer).
    // sum_layers --> all unique synthetic (layer, cellID) channels.
    for (const auto& [cell_id, stats] : layers[layer].channels) {
      (void)cell_id;
      const double avg = static_cast<double>(stats.total_hits) / static_cast<double>(n_events);
      const double rate = static_cast<double>(stats.total_hits) / (static_cast<double>(n_events) * kEventWindowSec);
      h_avg->Fill(stats.x, stats.y, avg);
      h_max->Fill(stats.x, stats.y, stats.max_hits_event);
      h_rate->Fill(stats.x, stats.y, rate);
      h_avg_r->Fill(stats.r, avg);
      h_nchan_r->Fill(stats.r, 1.0);
    }

    for (int bin = 1; bin <= h_avg_r->GetNbinsX(); ++bin) {
      const double nch = h_nchan_r->GetBinContent(bin);
      if (nch > 0.0) h_avg_r->SetBinContent(bin, h_avg_r->GetBinContent(bin) / nch);
    }

    // Output.
    auto* hist_dir = dir->mkdir("hists");
    hist_dir->cd();
    h_avg->Write();
    h_max->Write();
    h_rate->Write();
    h_avg_r->Write();
    h_nchan_r->Write();
    layers[layer].h_hits_evt->Write();

    draw_and_write(dir, h_avg, "c_avg", true, false);
    draw_and_write(dir, h_max, "c_max", true, false);
    draw_and_write(dir, h_rate, "c_rate", true, false);
    draw_and_write(dir, h_avg_r, "c_avg_r", false, false);
    draw_and_write(dir, layers[layer].h_hits_evt, "c_hits_evt", false, true);
  }

  std::cout << "\n";
  output.Close();
  return 0;
}
