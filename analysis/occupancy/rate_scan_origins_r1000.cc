/*

./build/rate_scan_origins_r1000 -i data/bkg_apr -o plots/occupancy/rate_scan_origins_r1000.root

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

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "classify_hit.h"
#include "decode_cell_id.h"
#include "utils.h"

#include <algorithm>
#include <array>
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
constexpr double MIP_2 = 7.0e-3;
constexpr std::array<double, 6> kCoefficients = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
const std::array<int, kCoefficients.size()> kColors = {kBlack, kBlue + 1, kGreen + 2, kOrange + 1, kRed + 1, kMagenta + 1};

struct OriginSpec {
  const char* label = "";
  int origin_index = -1;
};

constexpr std::array<OriginSpec, 2> kOrigins = {{
    {"DIS", 0},
    {"pBeamGas", 5},
}};

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
  double p99_rate_hz = 0.0;
};

struct EventChannel {
  double energy_gev = 0.0;
  std::array<double, 2> energy_by_origin{};
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

std::string percentile_label(double threshold_mip, double percentile_rate_hz, int percentile) {
  char buffer[128];
  std::snprintf(buffer,
                sizeof(buffer),
                "%s (p%d %.3g Hz)",
                threshold_label(threshold_mip).c_str(),
                percentile,
                percentile_rate_hz);
  return buffer;
}

void draw_rate_overlay(TDirectory* dir,
                       const std::array<ThresholdProducts, kCoefficients.size()>& products,
                       const std::string& origin_label,
                       int layer,
                       int percentile) {
  dir->cd();
  TCanvas canvas(("c_" + origin_label + "_layer" + std::to_string(layer) + "_rate_threshold_overlay_p" + std::to_string(percentile) + "_r1000").c_str(),
                 (origin_label + ", layer " + std::to_string(layer) + " rate threshold overlay, R < 1000 mm").c_str(),
                 1000,
                 800);
  canvas.SetLogx();
  canvas.SetLogy();

  TLegend legend(0.765, 0.65, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);

  double max_y = 0.0;

  // To set p95 line height later.
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
    hist->SetTitle((origin_label + ", layer " + std::to_string(layer) + ", R < 1000 mm;rate [Hz/channel];channels").c_str());
    hist->Draw(drew ? "hist same" : "hist");
    drew = true;
    const double percentile_rate_hz = percentile == 95 ? products[threshold_index].p95_rate_hz
                                                       : products[threshold_index].p99_rate_hz;
    legend.AddEntry(hist,
                    percentile_label(kCoefficients[threshold_index], percentile_rate_hz, percentile).c_str(),
                    "l");
  }

  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    const auto& product = products[threshold_index];
    const double percentile_rate_hz = percentile == 95 ? product.p95_rate_hz : product.p99_rate_hz;
    if (percentile_rate_hz <= 0.0) continue;

    auto* line = new TLine(percentile_rate_hz, kHistMinimum, percentile_rate_hz, 1.25 * max_y);
    line->SetLineColor(kColors[threshold_index]);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
  }

  legend.Draw();
  canvas.Write();
}

int tracked_origin_slot(int raw_origin_index) {
  if (raw_origin_index == kOrigins[0].origin_index) return 0;
  if (raw_origin_index == kOrigins[1].origin_index) return 1;
  return -1;
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and find ROOT files.
  const auto args = parse_args(argc, argv);
  const auto files = br::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  // Ensure output stuff.
  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // One rate histogram per origin, readout layer, and threshold.
  const auto rate_edges = br::log_edges(240, 1.0, 1e7);
  std::array<std::array<std::array<ThresholdProducts, kCoefficients.size()>, kNReadoutLayers>, kOrigins.size()> products{};

  for (std::size_t origin_slot = 0; origin_slot < kOrigins.size(); ++origin_slot) {
    for (int layer = 0; layer < kNReadoutLayers; ++layer) {
      for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
        const std::string name =
            "h_channel_rate_" + std::string(kOrigins[origin_slot].label) +
            "_layer" + std::to_string(layer) + "_" + threshold_tag(kCoefficients[threshold_index]) + "_r1000";
        const std::string title =
            "LFHCAL per-channel rate, " + std::string(kOrigins[origin_slot].label) +
            ", layer " + std::to_string(layer) +
            ", threshold " + threshold_label(kCoefficients[threshold_index]) +
            ", R < 1000 mm;rate [Hz/channel];channels";
        products[origin_slot][layer][threshold_index].h_rate = new TH1D(name.c_str(), title.c_str(), 240, rate_edges.data());
      }
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

      // For this event, store the summed channel signal before applying thresholds.
      std::array<std::unordered_map<br::LFHCALChannelID, EventChannel, br::LFHCALChannelIDHash>, kNReadoutLayers> event_channels;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop hits.
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        const auto channel = decoder.channel(cell_id);
        if (channel.rlayerz < 0 || channel.rlayerz >= kNReadoutLayers) continue;

        // Keep just the central channels whose transverse radius is below 1000 mm.
        const auto position = decoder.position(cell_id);
        const double radius_mm2 = position.x_mm * position.x_mm + position.y_mm * position.y_mm;
        if (radius_mm2 > kMaxRadiusMm2) continue;

        // Does one of two things:
        // 1. There are no preexisting event channels with this ID --> create new event channel.
        // 2. There is a preexisting event channel with this ID --> accumulate into that one.
        auto& event_channel = event_channels[channel.rlayerz][channel];
        event_channel.energy_gev += hit.getEnergy(); // Sum energy here to apply a threshold later.

        // Per hit, sum all contributions that fall into the tracked truth families.
        for (const auto& contribution : hit.getContributions()) {
          const double contribution_energy = contribution.getEnergy();
          if (contribution_energy <= 0.0) continue;

          const int raw_origin_index = br::origin_index(contribution.getParticle().getGeneratorStatus());
          const int origin_slot = tracked_origin_slot(raw_origin_index);
          if (origin_slot < 0) continue;
          event_channel.energy_by_origin[origin_slot] += contribution_energy;
        }
      }

      // Increment channels that pass each threshold for each tracked origin.
      for (int layer = 0; layer < kNReadoutLayers; ++layer) {
        for (const auto& [channel, event_channel] : event_channels[layer]) {
          for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

            // Apply the summed channel threshold.
            // *No thresholds applied at the contribution level: if a channel exceeds the threshold, all its 
            // contributions will make it*
            if (event_channel.energy_gev <= kCoefficients[threshold_index] * mip_energy_gev(layer)) continue;

            for (std::size_t origin_slot = 0; origin_slot < kOrigins.size(); ++origin_slot) {

              // Skip if this origin contributes no energy to the channel in this event.
              if (event_channel.energy_by_origin[origin_slot] <= 0.0) continue;
              ++products[origin_slot][layer][threshold_index].channel_passes[channel];
            }
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

  // Go through all tracked origins and readout layers.
  for (std::size_t origin_slot = 0; origin_slot < kOrigins.size(); ++origin_slot) {
    for (int layer = 0; layer < kNReadoutLayers; ++layer) {
      for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
        std::vector<double> rates_hz;

        // Loop through channels that pass this threshold for this origin.
        for (const auto& [channel, passes] : products[origin_slot][layer][threshold_index].channel_passes) {
          (void)channel;

          // Compute and fill channel rate into hist.
          const double rate_hz = static_cast<double>(passes) / total_time_sec;
          products[origin_slot][layer][threshold_index].h_rate->Fill(rate_hz);
          rates_hz.push_back(rate_hz);
        }

        if (!rates_hz.empty()) {
          std::sort(rates_hz.begin(), rates_hz.end());
          const std::size_t p95_index = static_cast<std::size_t>(0.95 * static_cast<double>(rates_hz.size() - 1));
          const std::size_t p99_index = static_cast<std::size_t>(0.99 * static_cast<double>(rates_hz.size() - 1));
          products[origin_slot][layer][threshold_index].p95_rate_hz = rates_hz[p95_index];
          products[origin_slot][layer][threshold_index].p99_rate_hz = rates_hz[p99_index];
        }
      }
    }
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  // Write hists and overlays under one top-level directory per truth family.
  for (std::size_t origin_slot = 0; origin_slot < kOrigins.size(); ++origin_slot) {
    auto* origin_dir = output.mkdir(kOrigins[origin_slot].label);
    auto* hist_dir = origin_dir->mkdir("hists");
    auto* p95_dir = origin_dir->mkdir("p95");
    auto* p99_dir = origin_dir->mkdir("p99");
    hist_dir->cd();

    for (int layer = 0; layer < kNReadoutLayers; ++layer) {
      for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
        products[origin_slot][layer][threshold_index].h_rate->Write();
      }
    }

    for (int layer = 0; layer < kNReadoutLayers; ++layer) {
      draw_rate_overlay(p95_dir, products[origin_slot][layer], kOrigins[origin_slot].label, layer, 95);
      draw_rate_overlay(p99_dir, products[origin_slot][layer], kOrigins[origin_slot].label, layer, 99);
    }
  }

  std::cout << "\n";
  output.Close();
  return 0;
}
