/*

./build/threshold_rate_overlays -i data/reco_bkg_apr -o plots/rate_correlations/threshold_rate_overlays.root

*/

#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "utils.h"

#include <array>
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

constexpr std::array<double, 5> kThresholdsGeV = {0.0, 0.001, 0.002, 0.003, 0.004};
const std::array<int, kThresholdsGeV.size()> kColors = {kBlack, kBlue + 1, kGreen + 2, kOrange + 1, kRed + 1};

struct Args {
  std::string input_dir;
  std::string output_file;
};

struct ThresholdProducts {
  TH1D* h_rate = nullptr;
  std::unordered_map<br::LFHCALChannelID, std::uint64_t, br::LFHCALChannelIDHash> channel_passes;
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

std::string threshold_label(double threshold_gev) {
  return std::to_string(static_cast<int>(threshold_gev * 1000.0 + 0.5)) + " MeV";
}

std::string threshold_tag(double threshold_gev) {
  return "thr" + std::to_string(static_cast<int>(threshold_gev * 1000.0 + 0.5)) + "MeV";
}

void draw_rate_overlay(TDirectory* dir,
                       const std::array<ThresholdProducts, kThresholdsGeV.size()>& products,
                       int layer) {
  dir->cd();
  TCanvas canvas(("c_layer" + std::to_string(layer) + "_rate_threshold_overlay").c_str(),
                 ("Layer " + std::to_string(layer) + " rate threshold overlay").c_str(),
                 1000,
                 800);
  canvas.SetLogx();
  canvas.SetLogy();

  TLegend legend(0.65, 0.65, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);

  bool drew = false;
  for (std::size_t threshold_index = 0; threshold_index < kThresholdsGeV.size(); ++threshold_index) {
    auto* hist = products[threshold_index].h_rate;
    hist->SetStats(false);
    hist->SetLineColor(kColors[threshold_index]);
    hist->SetLineWidth(2);
    hist->SetTitle(("Layer " + std::to_string(layer) + ";rate [Hz/channel];channels").c_str());
    hist->Draw(drew ? "hist same" : "hist");
    drew = true;
    legend.AddEntry(hist, threshold_label(kThresholdsGeV[threshold_index]).c_str(), "l");
  }

  legend.Draw();
  canvas.Write();
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  const auto args = parse_args(argc, argv);
  const auto files = br::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  const auto rate_edges = br::log_edges(240, 1.0, 1e7);
  std::array<std::array<ThresholdProducts, kThresholdsGeV.size()>, kNReadoutLayers> products{};

  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    for (std::size_t threshold_index = 0; threshold_index < kThresholdsGeV.size(); ++threshold_index) {
      const std::string name = "h_channel_rate_layer" + std::to_string(layer) + "_" + threshold_tag(kThresholdsGeV[threshold_index]);
      const std::string title = "LFHCAL per-channel rate, layer " + std::to_string(layer) +
          ", threshold " + threshold_label(kThresholdsGeV[threshold_index]) +
          ";rate [Hz/channel];channels";
      products[layer][threshold_index].h_rate = new TH1D(name.c_str(), title.c_str(), 240, rate_edges.data());
    }
  }

  const br::LFHCALCellIDDecoder decoder;
  std::uint64_t n_events = 0;
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
      if (!br::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      std::array<std::unordered_map<br::LFHCALChannelID, double, br::LFHCALChannelIDHash>, kNReadoutLayers> channel_energy_by_layer;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop hits and sum deposited energy per channel in this event.
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        const auto channel = decoder.channel(cell_id);
        if (channel.rlayerz < 0 || channel.rlayerz >= kNReadoutLayers) continue;

        channel_energy_by_layer[channel.rlayerz][channel] += hit.getEnergy();
      }

      // Increment channels that pass each threshold for this event.
      for (int layer = 0; layer < kNReadoutLayers; ++layer) {
        for (const auto& [channel, energy] : channel_energy_by_layer[layer]) {
          for (std::size_t threshold_index = 0; threshold_index < kThresholdsGeV.size(); ++threshold_index) {
            if (energy <= kThresholdsGeV[threshold_index]) continue;
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

  // For each channel, calculate: total passes / total time
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    for (std::size_t threshold_index = 0; threshold_index < kThresholdsGeV.size(); ++threshold_index) {
      for (const auto& [channel, passes] : products[layer][threshold_index].channel_passes) {
        (void)channel;
        products[layer][threshold_index].h_rate->Fill(static_cast<double>(passes) / total_time_sec);
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
    for (std::size_t threshold_index = 0; threshold_index < kThresholdsGeV.size(); ++threshold_index) {
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
