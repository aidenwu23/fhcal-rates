/*

./build/rate_vs_mip_r1000 -i data/bkg_apr -o plots/occupancy/rate_vs_mip_r1000.root

*/

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kHitCollection = "LFHCALHits";
constexpr int kNReadoutLayers = 7;
constexpr double kEventWindowSec = 2e-6;
constexpr double kMaxRadiusMm = 1000.0;
constexpr double kMaxRadiusMm2 = kMaxRadiusMm * kMaxRadiusMm;

// The first two readout layers use one MIP scale, and the later layers use another.
constexpr double MIP_1 = 3.5e-3;
constexpr double MIP_2 = 7.0e-3;
constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
constexpr std::array<double, 2> kPercentiles = {0.95, 0.99};
constexpr std::array<int, kNReadoutLayers> kLayerColors = {
    kBlack, kBlue + 1, kGreen + 2, kOrange + 1, kRed + 1, kMagenta + 1, kCyan + 1};

struct Args {
  std::string input_dir;
  std::string output_file;
};

struct ThresholdProducts {
  // channel_passes[channel]: number of events where this channel cleared one coefficient choice.
  std::unordered_map<br::LFHCALChannelID, std::uint64_t, br::LFHCALChannelIDHash> channel_passes;
};

// ----------------------------------------------------------------------------------
// CLI stuff.
// ----------------------------------------------------------------------------------
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
// Label helpers and drawing.
// ----------------------------------------------------------------------------------
// Convert a layer index into the MIP energy scale used for that layer family.
double mip_energy_gev(int layer) {
  return layer < 2 ? MIP_1 : MIP_2;
}

// Turn a percentile request into the corresponding data index.
std::size_t percentile_index(std::size_t n_values, double percentile) {
  if (n_values == 0) return 0;
  return static_cast<std::size_t>(percentile * static_cast<double>(n_values - 1));
}

std::string percentile_tag(double percentile) {
  return "p" + std::to_string(static_cast<int>(percentile * 100.0 + 0.5));
}

std::string percentile_title(double percentile) {
  return percentile_tag(percentile) + " rate vs MIP coefficient, R < 1000 mm;MIP coefficient;rate [Hz/channel]";
}

// Write one overlay plot for one percentile choice.
// Each curve is one readout layer, and the x axis scans the MIP coefficient.
void draw_overlay(TFile& output,
                  double percentile,
                  const std::array<std::array<double, kCoefficients.size()>, kNReadoutLayers>& values) {
  // Keep p95 and p99 products in separate top-level folders.
  auto* dir = output.mkdir(percentile_tag(percentile).c_str());
  dir->cd();

  TCanvas canvas(("c_" + percentile_tag(percentile) + "_rate_vs_mip").c_str(),
                 percentile_title(percentile).c_str(),
                 1000,
                 800);
  canvas.SetLogy();

  // Find the tallest point first so every layer fits on the same frame.
  double max_y = 0.0;
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    for (double value : values[layer]) {
      max_y = std::max(max_y, value);
    }
  }
  if (max_y <= 0.0) max_y = 1.0;

  auto* frame = canvas.DrawFrame(kCoefficients.front(), 1.0, kCoefficients.back(), 1.25 * max_y);
  frame->SetTitle(percentile_title(percentile).c_str());
  frame->GetXaxis()->SetTitle("MIP coefficient");
  frame->GetYaxis()->SetTitle("rate [Hz/channel]");
  frame->SetStats(false);

  TLegend legend(0.65, 0.60, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);

  // Build one graph per layer and draw them onto the shared frame.
  std::array<TGraph*, kNReadoutLayers> graphs{};
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    graphs[layer] = new TGraph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), values[layer].data());
    graphs[layer]->SetName(("g_" + percentile_tag(percentile) + "_layer" + std::to_string(layer)).c_str());
    graphs[layer]->SetTitle(("Layer " + std::to_string(layer)).c_str());
    graphs[layer]->SetLineColor(kLayerColors[layer]);
    graphs[layer]->SetMarkerColor(kLayerColors[layer]);
    graphs[layer]->SetMarkerStyle(20 + layer);
    graphs[layer]->SetLineWidth(2);
    graphs[layer]->Draw("LP SAME");
    legend.AddEntry(graphs[layer], ("Layer " + std::to_string(layer)).c_str(), "lp");
  }

  legend.Draw();
  canvas.Write();
  for (auto* graph : graphs) graph->Write();
  output.cd();
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse inputs and make sure there is data to process.
  const auto args = parse_args(argc, argv);
  const auto files = br::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // products[layer][coefficient]: persistent event-pass counts for that layer and threshold choice.
  std::array<std::array<ThresholdProducts, kCoefficients.size()>, kNReadoutLayers> products{};
  const br::LFHCALCellIDDecoder decoder;
  std::uint64_t n_events = 0;
  br::FileProgress progress(files.size(), std::cerr);

  // Loop over all files and accumulate channel pass counts event by event.
  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!br::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      // For this event, store the summed channel energy separately in each readout layer.
      std::array<std::unordered_map<br::LFHCALChannelID, double, br::LFHCALChannelIDHash>, kNReadoutLayers> channel_energy_by_layer;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop hits and fold them into the event-level channel sums.
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        const auto channel = decoder.channel(cell_id);
        if (channel.rlayerz < 0 || channel.rlayerz >= kNReadoutLayers) continue;

        const auto position = decoder.position(cell_id);
        const double radius_mm2 = position.x_mm * position.x_mm + position.y_mm * position.y_mm;
        if (radius_mm2 > kMaxRadiusMm2) continue;

        channel_energy_by_layer[channel.rlayerz][channel] += hit.getEnergy();
      }

      // After the event is summed, count which channels passed each coefficient threshold.
      for (int layer = 0; layer < kNReadoutLayers; ++layer) {
        for (const auto& [channel, energy] : channel_energy_by_layer[layer]) {
          for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
            if (energy <= kCoefficients[threshold_index] * mip_energy_gev(layer)) continue;
            ++products[layer][threshold_index].channel_passes[channel];
          }
        }
      }
    }
  }
  std::cerr << "\n";

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  // Convert event counts into rates, then extract the requested percentiles for each layer.
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<std::array<std::array<double, kCoefficients.size()>, kNReadoutLayers>, kPercentiles.size()> percentile_values{};

  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
      std::vector<double> rates_hz;
      // Collect this layer-and-threshold rate from every channel that ever passed.
      for (const auto& [channel, passes] : products[layer][threshold_index].channel_passes) {
        (void)channel;
        rates_hz.push_back(static_cast<double>(passes) / total_time_sec);
      }
      if (rates_hz.empty()) continue;

      // Sort once so percentile lookup becomes a simple indexed read.
      std::sort(rates_hz.begin(), rates_hz.end());
      for (std::size_t percentile_index_slot = 0; percentile_index_slot < kPercentiles.size(); ++percentile_index_slot) {
        const std::size_t index = percentile_index(rates_hz.size(), kPercentiles[percentile_index_slot]);
        percentile_values[percentile_index_slot][layer][threshold_index] = rates_hz[index];
      }
    }
  }

  // Write one overlay under p95 and one under p99.
  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  for (std::size_t percentile_index_slot = 0; percentile_index_slot < kPercentiles.size(); ++percentile_index_slot) {
    draw_overlay(output, kPercentiles[percentile_index_slot], percentile_values[percentile_index_slot]);
  }

  output.Close();
  return 0;
}
