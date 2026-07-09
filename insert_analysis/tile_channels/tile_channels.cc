/*

./build/insert_tile_channels -i data/bkg_apr -o insert_plots/tile_channels.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "insert_to_lfhcal.h"
#include "utils.h"

#include <algorithm>
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

// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
constexpr const char* kInsertCollectionMatch = "HcalEndcapPInsert";
constexpr double kEventWindowSec = 2e-6;
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;
constexpr double kTileMipGeV = 4e-4;
constexpr std::array<double, 16> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
constexpr std::array<double, 2> kPercentiles = {0.95, 0.99};
const std::array<int, 2> kColors = {kBlue + 1, kRed + 1};
const std::array<int, 2> kMarkers = {20, 21};

struct Args {
  std::string input_path;
  std::string output_file;
};

struct ThresholdSum {
  std::unordered_map<rates::VirtualLFHCALChannelID, std::uint64_t, rates::VirtualLFHCALChannelIDHash> pass_counts;
  std::unordered_map<rates::VirtualLFHCALChannelID, double, rates::VirtualLFHCALChannelIDHash> payload_bits;
};

using ChannelEnergyMap = std::unordered_map<rates::VirtualLFHCALChannelID, double, rates::VirtualLFHCALChannelIDHash>;

// ----------------------------------------------------------------------------------
// CLI handling.
// ----------------------------------------------------------------------------------
void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT.root|DIR -o OUTPUT.root\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_path = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      args.output_file = argv[++i];
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input_path.empty() || args.output_file.empty()) {
    usage(argv[0]);
    std::exit(1);
  }

  return args;
}

// ----------------------------------------------------------------------------------
// Helpers and plotters.
// ----------------------------------------------------------------------------------
std::string find_insert_collection(const podio::Frame& frame) {
  for (const auto& name : frame.getAvailableCollections()) {
    if (name.find(kInsertCollectionMatch) != std::string::npos &&
        name.find("Contributions") == std::string::npos) {
      return name;
    }
  }
  return "";
}

std::size_t percentile_index(std::size_t n_values, double percentile) {
  if (n_values == 0) return 0;
  return static_cast<std::size_t>(percentile * static_cast<double>(n_values - 1));
}

std::string percentile_tag(double percentile) {
  return "p" + std::to_string(static_cast<int>(percentile * 100.0 + 0.5));
}

void draw_overlay(TDirectory* dir,
                  const char* canvas_name,
                  const char* title,
                  const char* graph_prefix,
                  const char* y_title,
                  double y_min,
                  double y_max,
                  const std::array<std::array<double, 16>, 2>& percentile_values) {
  dir->cd();

  TCanvas canvas(canvas_name, title, 1000, 800);
  canvas.SetGrid();

  auto* frame = canvas.DrawFrame(kCoefficients.front(), y_min, kCoefficients.back(), y_max);
  frame->SetTitle(title);
  frame->GetXaxis()->SetTitle("MIP coefficient");
  frame->GetYaxis()->SetTitle(y_title);
  frame->SetStats(false);

  TLegend legend(0.65, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.04);

  std::array<TGraph, 2> graphs = {
      TGraph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), percentile_values[0].data()),
      TGraph(static_cast<int>(kCoefficients.size()), kCoefficients.data(), percentile_values[1].data())};

  for (int i = static_cast<int>(kPercentiles.size()) - 1; i >= 0; --i) {
    graphs[i].SetName((std::string(graph_prefix) + "_" + percentile_tag(kPercentiles[i])).c_str());
    graphs[i].SetLineColor(kColors[i]);
    graphs[i].SetMarkerColor(kColors[i]);
    graphs[i].SetMarkerStyle(kMarkers[i]);
    graphs[i].SetLineWidth(2);
    graphs[i].Draw("LP SAME");
    legend.AddEntry(&graphs[i], percentile_tag(kPercentiles[i]).c_str(), "lp");
  }

  legend.Draw();
  canvas.Write();
  for (auto& graph : graphs) graph.Write();
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main.
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  const auto args = parse_args(argc, argv);
  std::vector<fs::path> input_files;
  const fs::path input_path(args.input_path);
  if (fs::is_regular_file(input_path) && input_path.extension() == ".root") {
    input_files.push_back(input_path);
  } else {
    input_files = rates::find_root_files(args.input_path);
  }
  if (input_files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_path << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  const rates::HcalEndcapPInsertCellIDDecoder decoder;
  const rates::InsertToLFHCALMapper mapper;

  std::array<ThresholdSum, 16> threshold_sums;
  std::string insert_collection_name;
  std::uint64_t n_events = 0;
  rates::FileProgress progress(input_files.size(), std::cerr);

  // Loop over all input filels.
  for (const auto& path : input_files) {
    progress.tick();

    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Per file, loop over all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto event_data = reader.readEntry("events", event_index);
      if (!event_data) continue;

      podio::Frame frame(std::move(event_data));
      if (insert_collection_name.empty()) {
        insert_collection_name = find_insert_collection(frame);
      }
      if (insert_collection_name.empty() || !rates::has_collection(frame, insert_collection_name)) continue;
      ++n_events;

      ChannelEnergyMap channel_energies;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(insert_collection_name);

      // Per event, loop over all hits.
      for (const auto& hit : hits) {

        // Decode hit onto cell.
        const auto cell = decoder.cell(static_cast<std::uint64_t>(hit.getCellID()));
        if (cell.layer < 1 || cell.layer > 60) continue;

        // Map cell into a virtual tile based on position and treat it as a channel.
        const auto position = hit.getPosition();
        auto channel = mapper.channel(cell.layer, position.x, position.y);
        channel.layer = cell.layer;
        channel_energies[channel] += hit.getEnergy();
      }

      // After accumulating all hits in this event, loop over all channels.
      for (const auto& [channel, energy] : channel_energies) {

        // Per channel, loop over all thresholds.
        for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

          // Skip if channel energy doeosn't exceed threshold.
          if (energy <= kCoefficients[threshold_index] * kTileMipGeV) continue;

          // Otherwise increment corresponding stats.
          ++threshold_sums[threshold_index].pass_counts[channel];
          threshold_sums[threshold_index].payload_bits[channel] += (kOverheadBits + kBitsPerHit) * kSamplesPerEvent;
        }
      }
    }
  }
  std::cerr << "\n";

  if (insert_collection_name.empty()) {
    std::cerr << "Failed to find an insert hit collection\n";
    return 1;
  }
  if (n_events == 0) {
    std::cerr << "No events with insert hits found\n";
    return 1;
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::array<std::array<double, 16>, 2> channel_percentile_values{};
  std::array<std::array<double, 16>, 2> data_percentile_values{};

  // Loop over all thresholds.
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    std::vector<double> channel_rates_hz;

    // Per threshold, loop over all passed channels.
    for (const auto& [channel, pass_count] : threshold_sums[threshold_index].pass_counts) {
      (void)channel;

      // Calc and store stats for this channel.
      channel_rates_hz.push_back(static_cast<double>(pass_count) / total_time_sec);
    }
    if (!channel_rates_hz.empty()) {

      // Sort all channel stats.
      std::sort(channel_rates_hz.begin(), channel_rates_hz.end());

      // Grab the requested percentiles.
      for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
        channel_percentile_values[percentile_slot][threshold_index] =
            channel_rates_hz[percentile_index(channel_rates_hz.size(), kPercentiles[percentile_slot])];
      }
    }

    // Do the same but with Gb/s data rate.
    std::vector<double> data_rates_gbps;
    for (const auto& [channel, payload_bits] : threshold_sums[threshold_index].payload_bits) {
      (void)channel;
      data_rates_gbps.push_back((payload_bits / total_time_sec) / 1.0e9);
    }
    if (!data_rates_gbps.empty()) {
      std::sort(data_rates_gbps.begin(), data_rates_gbps.end());
      for (std::size_t percentile_slot = 0; percentile_slot < kPercentiles.size(); ++percentile_slot) {
        data_percentile_values[percentile_slot][threshold_index] =
            data_rates_gbps[percentile_index(data_rates_gbps.size(), kPercentiles[percentile_slot])];
      }
    }
  }

  auto* channel_dir = output.mkdir("channel_rate_vs_mip");
  draw_overlay(channel_dir,
               "c_channel_rate_vs_mip",
               "Tile channel tail rates vs MIP coefficient;MIP coefficient;rate [Hz]",
               "g_channel_rate_vs_mip",
               "rate [Hz]",
               0.0,
               7.0e4,
               channel_percentile_values);

  auto* data_dir = output.mkdir("data_rate_vs_mip");
  draw_overlay(data_dir,
               "c_data_rate_vs_mip",
               "Tile channel tail data rates vs MIP coefficient;MIP coefficient;data rate [Gb/s]",
               "g_data_rate_vs_mip",
               "data rate [Gb/s]",
               0.0,
               4.5e-2,
               data_percentile_values);

  output.Close();
  return 0;
}
