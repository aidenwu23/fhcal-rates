/*

./build/occupancy -i data/reco_bkg_feb -o plots/occupancy/occupancy.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4eic/CalorimeterHitCollection.h>

#include "decode_channel.h"
#include "smooth_hists.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <unordered_set>
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
constexpr double kCellSizeMM = 51.2 ;                    // Display bin size in x and y.
constexpr double kXYExtentMM = 2700.0;                  // +- 2700 mm for x and y ranges.

// Stats per channel (cells sharing transverse location in a readout layer).
struct ChannelStats {
  std::unordered_set<std::uint64_t> raw_cell_ids; // Keep track of which cells belong in this channel.

  // Accumulate the raw hit positions from all cells in this channel.
  double x_sum = 0.0;
  double y_sum = 0.0;

  // Total number of positions added.
  std::uint64_t n_pos = 0;
  std::uint64_t total_hits = 0;

  int max_hits_event = 0; // Keep track of which event had the most hits. 

  // Divide sum of positions by total positions added to get the average hit positions for the whole channel.
  double x() const { return n_pos > 0 ? x_sum / static_cast<double>(n_pos) : 0.0; }
  double y() const { return n_pos > 0 ? y_sum / static_cast<double>(n_pos) : 0.0; }
  double r() const { return std::hypot(x(), y()); } // radial distance from center.
};

// Hold one layer's data.
struct LayerAccum {

  // For each unique channelID, it stores the statistics and also the ID hash.
  std::unordered_map<br::LFHCALChannelID, ChannelStats, br::LFHCALChannelIDHash> channels;
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
  std::vector<LayerAccum> layers(kNLayers + 1); // +1 to account for the summed layer stats.
  for (int layer = 0; layer <= kNLayers; ++layer) {
    layers[layer].h_hits_evt = new TH1D(
        "h_hits_evt",
        (layer == kAllLayersIndex ? std::string("LFHCAL summed layers;hits/event;Events")
          : std::string("LFHCAL layer ") + std::to_string(layer) + ";hits/event;Events").c_str(),
        200, 0, 200);
  }

  std::uint64_t n_events = 0;

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

      // Grab LFHCALRecHits when possible.
      podio::Frame frame(std::move(data));
      if (!br::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      // Per-event channel counters.
      // event_counts[layer][channel]: For this event, the number of hits that landed in [layer]'s [channel].
      std::vector<std::unordered_map<br::LFHCALChannelID, int, br::LFHCALChannelIDHash>> event_counts(kNLayers + 1);
      std::vector<int> layer_totals(kNLayers + 1, 0);

      // Grab and loop over hits.
      const auto& hits = frame.get<edm4eic::CalorimeterHitCollection>(kHitCollection);
      for (const auto& hit : hits) {

        // Apply hit threshold.
        if (hit.getEnergy() <= threshold_geV) continue;

        // Get layer and cellID.
        const int layer = hit.getLayer();
        if (layer < 0 || layer >= kNLayers) continue;

        // Decode channel ID from cell ID.
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        const auto channel_id = decoder.channel(cell_id);
        const auto pos = hit.getPosition();

        // Look up the global statistics object for this channel.
        // If this channel has never been seen before, the map creates a fresh ChannelStats entry, otherwise,
        // it conveniently grabs the preexisting one.
        auto& stats = layers[layer].channels[channel_id];
        stats.raw_cell_ids.insert(cell_id);
        stats.x_sum += pos.x;
        stats.y_sum += pos.y;
        ++stats.n_pos;

        // Count one more hit in this channel for the current event.
        ++event_counts[layer][channel_id];
        ++layer_totals[layer];

        // Also update the "all layers together" bucket.
        auto& merged_stats = layers[kAllLayersIndex].channels[channel_id];
        merged_stats.raw_cell_ids.insert(cell_id);
        merged_stats.x_sum += pos.x;
        merged_stats.y_sum += pos.y;
        ++merged_stats.n_pos;
        ++event_counts[kAllLayersIndex][channel_id];
        ++layer_totals[kAllLayersIndex];
      }

      // End of this event:
      // move the temporary per-event counts into the long-lived channel statistics.
      for (int layer = 0; layer <= kNLayers; ++layer) {
        layers[layer].h_hits_evt->Fill(layer_totals[layer]);
        for (const auto& [channel_id, count] : event_counts[layer]) {
          auto& stats = layers[layer].channels[channel_id];

          // Add this event's contribution to the all-events total.
          stats.total_hits += static_cast<std::uint64_t>(count);

          // Keep the largest number of hits this channel ever saw in any single event.
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
    for (const auto& [channel_id, stats] : layers[layer].channels) {
      (void)channel_id;
      rmax = std::max(rmax, stats.r());
    }
    const double xmin = -kXYExtentMM;
    const double xmax = kXYExtentMM;
    const double ymin = -kXYExtentMM;
    const double ymax = kXYExtentMM;
    const int xbins = static_cast<int>(std::ceil((xmax - xmin) / kCellSizeMM));
    const int ybins = static_cast<int>(std::ceil((ymax - ymin) / kCellSizeMM));

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
    for (const auto& [channel_id, stats] : layers[layer].channels) {
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
