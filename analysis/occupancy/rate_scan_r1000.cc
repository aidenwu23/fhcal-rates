/*

./build/rate_scan_r1000 -i data/reco_bkg_apr -o plots/occupancy/rate_scan_r1000.root

*/

#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
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
constexpr double kHistMinimum = 0.8;

constexpr double MIP_1 = 3.5e-3;
constexpr double MIP_2 = 7.25e-3;
constexpr std::array<double, 11> kCoefficients = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
const std::array<int, kCoefficients.size()> kColors = {kBlack, kBlue + 1, kGreen + 2, kOrange + 1, kRed + 1};

// ----------------------------------------------------------------------------------
// CLI and per-threshold products.
// ----------------------------------------------------------------------------------
struct Args {
  std::string input_dir;
  std::string output_file;
};

struct ThresholdProducts {
  TH1D* h_rate = nullptr;
  std::unordered_map<br::LFHCALChannelID, std::uint64_t, br::LFHCALChannelIDHash> channel_passes;
  double p95_rate_hz = 0.0;
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
// Label helpers and drawing.
// ----------------------------------------------------------------------------------
double mip_energy_gev(int layer) {
  return layer < 2 ? MIP_1 : MIP_2;
}

std::string threshold_label(double threshold_mip) {
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "%.2f MIP", threshold_mip);
  return buffer;
}

std::string threshold_tag(double threshold_mip) {
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "thr%03dMIP", static_cast<int>(threshold_mip * 100.0 + 0.5));
  return buffer;
}

std::string percentile_label(double threshold_mip, double p95_rate_hz) {
  char buffer[128];
  std::snprintf(buffer,
                sizeof(buffer),
                "%s (p95 %.3g Hz)",
                threshold_label(threshold_mip).c_str(),
                p95_rate_hz);
  return buffer;
}

void draw_rate_overlay(TDirectory* dir,
                       const std::array<ThresholdProducts, kCoefficients.size()>& products,
                       int layer) {
  dir->cd();
  TCanvas canvas(("c_layer" + std::to_string(layer) + "_rate_threshold_overlay_r1000").c_str(),
                 ("Layer " + std::to_string(layer) + " rate threshold overlay, R < 1000 mm").c_str(),
                 1000,
                 800);
  canvas.SetLogx();
  canvas.SetLogy();

  TLegend legend(0.52, 0.60, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);

  double max_y = 0.0;
  for (const auto& product : products) {
    max_y = std::max(max_y, product.h_rate->GetMaximum());
  }
  if (max_y <= 0.0) max_y = 1.0;

  bool drew = false;
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    auto* hist = products[threshold_index].h_rate;
    hist->SetStats(false);
    hist->SetLineColor(kColors[threshold_index]);
    hist->SetLineWidth(2);
    hist->SetMinimum(kHistMinimum);
    hist->SetMaximum(1.25 * max_y);
    hist->SetTitle(("Layer " + std::to_string(layer) + ", R < 1000 mm;rate [Hz/channel];channels").c_str());
    hist->Draw(drew ? "hist same" : "hist");
    drew = true;
    legend.AddEntry(hist,
                    percentile_label(kCoefficients[threshold_index], products[threshold_index].p95_rate_hz).c_str(),
                    "l");
  }

  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    const auto& product = products[threshold_index];
    if (product.p95_rate_hz <= 0.0) continue;

    auto* line = new TLine(product.p95_rate_hz, kHistMinimum, product.p95_rate_hz, 1.25 * max_y);
    line->SetLineColor(kColors[threshold_index]);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
  }

  legend.Draw();
  canvas.Write();
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args.
  const auto args = parse_args(argc, argv);
  const auto files = br::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // One rate histogram per readout layer and threshold.
  const auto rate_edges = br::log_edges(240, 1.0, 1e7);
  std::array<std::array<ThresholdProducts, kCoefficients.size()>, kNReadoutLayers> products{};

  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
      const std::string name = "h_channel_rate_layer" + std::to_string(layer) + "_" + threshold_tag(kCoefficients[threshold_index]) + "_r1000";
      const std::string title = "LFHCAL per-channel rate, layer " + std::to_string(layer) +
          ", threshold " + threshold_label(kCoefficients[threshold_index]) +
          ", R < 1000 mm;rate [Hz/channel];channels";
      products[layer][threshold_index].h_rate = new TH1D(name.c_str(), title.c_str(), 240, rate_edges.data());
    }
  }

  const br::LFHCALCellIDDecoder decoder;
  std::uint64_t n_events = 0;
  br::FileProgress progress(files.size(), std::cerr);

  // Loop input files.
  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Loop events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!br::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      // For this event, store the summed energy seen by each channel in each readout layer.
      std::array<std::unordered_map<br::LFHCALChannelID, double, br::LFHCALChannelIDHash>, kNReadoutLayers> channel_energy_by_layer;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop hits and sum deposited energy per channel in this event.
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        const auto channel = decoder.channel(cell_id);
        if (channel.rlayerz < 0 || channel.rlayerz >= kNReadoutLayers) continue;

        // Keep just the central channels whose transverse radius is below 1000 mm.
        const auto position = decoder.position(cell_id);
        const double radius_mm2 = position.x_mm * position.x_mm + position.y_mm * position.y_mm;
        if (radius_mm2 > kMaxRadiusMm2) continue;

        channel_energy_by_layer[channel.rlayerz][channel] += hit.getEnergy();
      }

      // Increment channels that pass each threshold for this event.
      for (int layer = 0; layer < kNReadoutLayers; ++layer) {
        for (const auto& [channel, energy] : channel_energy_by_layer[layer]) {
          for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
            if (energy <= kCoefficients[threshold_index] * mip_energy_gev(layer)) continue;

            // Count one threshold-passing event for this channel.
            ++products[layer][threshold_index].channel_passes[channel];
          }
        }
      }
    }
  }

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  // Go thru all rlayerz's.
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
      std::vector<double> rates_hz;

      // Loop thru channels that pass the threshold for each threshold.
      for (const auto& [channel, passes] : products[layer][threshold_index].channel_passes) {
        (void)channel;

        // Compute rate.
        const double rate_hz = static_cast<double>(passes) / total_time_sec;
        products[layer][threshold_index].h_rate->Fill(rate_hz);
        rates_hz.push_back(rate_hz);
      }

      if (!rates_hz.empty()) {
        std::sort(rates_hz.begin(), rates_hz.end());
        
        // Take the 95th quantile from all rates.
        const std::size_t index = static_cast<std::size_t>(0.95 * static_cast<double>(rates_hz.size() - 1));
        products[layer][threshold_index].p95_rate_hz = rates_hz[index];
      }
    }
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  auto* hist_dir = output.mkdir("hists");
  hist_dir->cd();
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
      products[layer][threshold_index].h_rate->Write();
    }
  }

  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    draw_rate_overlay(&output, products[layer], layer);
  }

  std::cout << "\n";
  output.Close();
  return 0;
}
