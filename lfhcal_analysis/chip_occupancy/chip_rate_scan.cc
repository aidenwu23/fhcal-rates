/*

./build/chip_rate_scan -i data/bkg_apr -o plots/chip_occupancy/chip_rate_scan.root

*/

#include <TCanvas.h>
#include <TColor.h>
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
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace {
// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
constexpr const char* kHitCollection = "LFHCALHits";
constexpr int kNReadoutLayers = 7;
constexpr double kEventWindowSec = 2e-6;
constexpr double kHistMinimum = 0.8;
constexpr std::array<double, 3> kCoefficients = {0.1, 0.5, 1.5};
const std::array<int, kCoefficients.size()> kColors = {kBlue + 1, kBlack, kRed + 1};

struct ThresholdProducts {
  TH1D* h_rate = nullptr;
  std::unordered_map<rates::LFHCALChipID, std::uint64_t, rates::LFHCALChipIDHash> chip_passes;
  double p95_rate_hz = 0.0;
  double p99_rate_hz = 0.0;
  double max_rate_hz = 0.0;
};

struct EventChip {
  std::unordered_map<rates::LFHCALChannelID, double, rates::LFHCALChannelIDHash> channel_energy;
};

// ----------------------------------------------------------------------------------
// CLI
// ----------------------------------------------------------------------------------
void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT_DIR -o OUTPUT.root\n";
}

struct Args {
  std::string input_dir;
  std::string output_file;
};

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
// Plot helper(s)
// ----------------------------------------------------------------------------------
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

std::string rate_label(const char* name, double rate_hz) {
  char buffer[128];
  std::snprintf(buffer, sizeof(buffer), "%s = %.3g Hz", name, rate_hz);
  return buffer;
}

std::size_t percentile_index(std::size_t n_values, double percentile) {
  if (n_values == 0) return 0;
  return static_cast<std::size_t>(percentile * static_cast<double>(n_values - 1));
}

void draw_single(TFile& output, const ThresholdProducts& product, double threshold_mip) {
  output.cd();

  TCanvas canvas(("c_chip_rate_" + threshold_tag(threshold_mip)).c_str(),
                 ("Chip rate histogram at " + threshold_label(threshold_mip) + " threshold;rate [Hz];chips").c_str(),
                 1000,
                 800);
  canvas.SetLogx();
  canvas.SetLogy();

  auto* hist = product.h_rate;
  hist->SetStats(false);
  hist->SetLineWidth(2);
  hist->SetMinimum(kHistMinimum);
  hist->SetTitle(("Chip rate histogram at " + threshold_label(threshold_mip) + " threshold;rate [Hz];chips").c_str());
  hist->Draw("hist");

  const double max_y = hist->GetMaximum() > 0.0 ? 1.25 * hist->GetMaximum() : 1.0;
  auto draw_marker = [&](double x_value, int color) {
    if (x_value <= 0.0) return;
    auto* line = new TLine(x_value, kHistMinimum, x_value, max_y);
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
  };

  draw_marker(product.p95_rate_hz, kBlue + 1);
  draw_marker(product.p99_rate_hz, kRed + 1);
  draw_marker(product.max_rate_hz, kGreen + 2);

  TLegend legend(0.62, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  TLine p95_line;
  p95_line.SetLineColor(kBlue + 1);
  p95_line.SetLineStyle(2);
  p95_line.SetLineWidth(2);
  legend.AddEntry(&p95_line, rate_label("p95", product.p95_rate_hz).c_str(), "l");
  TLine p99_line;
  p99_line.SetLineColor(kRed + 1);
  p99_line.SetLineStyle(2);
  p99_line.SetLineWidth(2);
  legend.AddEntry(&p99_line, rate_label("p99", product.p99_rate_hz).c_str(), "l");
  TLine max_line;
  max_line.SetLineColor(kGreen + 2);
  max_line.SetLineStyle(2);
  max_line.SetLineWidth(2);
  legend.AddEntry(&max_line, rate_label("max", product.max_rate_hz).c_str(), "l");
  legend.Draw();
  canvas.Write();
}

void draw_overlay(TFile& output, const std::array<ThresholdProducts, kCoefficients.size()>& products) {
  output.cd();

  TCanvas canvas("c_chip_rate_threshold_overlay",
                 "Chip rate distributions;rate [Hz/chip];chips",
                 1000,
                 800);
  canvas.SetLogx();
  canvas.SetLogy();

  TLegend legend(0.68, 0.72, 0.88, 0.88);
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
    hist->SetTitle("Chip rate distributions;rate [Hz/chip];chips");
    hist->SetLineColor(kColors[threshold_index]);
    hist->SetLineWidth(2);
    hist->SetMinimum(kHistMinimum);
    hist->SetMaximum(1.25 * max_y);
    hist->Draw(drew ? "hist same" : "hist");
    drew = true;
    legend.AddEntry(hist, threshold_label(kCoefficients[threshold_index]).c_str(), "l");
  }

  legend.Draw();
  canvas.Write();
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Grab args and files.
  const auto args = parse_args(argc, argv);
  const auto files = rates::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  const auto rate_edges = rates::log_edges(240, 1.0, 1e7);
  std::array<ThresholdProducts, kCoefficients.size()> products{};
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    const std::string name = "h_chip_rate_" + threshold_tag(kCoefficients[threshold_index]);
    const std::string title = "LFHCAL per-chip rate, threshold " + threshold_label(kCoefficients[threshold_index]) +
        ";rate [Hz/chip];chips";
    products[threshold_index].h_rate = new TH1D(name.c_str(), title.c_str(), 240, rate_edges.data());
  }

  const rates::LFHCALCellIDDecoder decoder;
  std::uint64_t n_events = 0;
  rates::FileProgress progress(files.size(), std::cerr);

  // Loop all files.
  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // For each file, loop all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!rates::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      std::unordered_map<rates::LFHCALChipID, EventChip, rates::LFHCALChipIDHash> event_chips;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // For each event, loop all hits.
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        // Increment the corresponding channel's energy for the corresponding readout chip.
        const auto channel = decoder.channel(cell_id);
        auto& event_chip = event_chips[decoder.decode_chip(cell_id)];

        event_chip.channel_energy[channel] += hit.getEnergy();
      }

      std::array<std::unordered_set<rates::LFHCALChipID, rates::LFHCALChipIDHash>, kCoefficients.size()> fired_chips;

      // After processing all hits into corresponding channels and chips, loop over all chips.
      for (const auto& [chip, event_chip] : event_chips) {

        // Per chip, loop over all channels.
        for (const auto& [channel, energy_gev] : event_chip.channel_energy) {
          const double mip_gev = rates::mip_energy_gev(channel.rlayerz);

          // Per channel, loop over all thresholds.
          for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

            // If channel energy meets threshold, count chip as fired.
            if (energy_gev <= kCoefficients[threshold_index] * mip_gev) continue;
            fired_chips[threshold_index].insert(chip);
          }
        }
      }

      for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
        for (const auto& chip : fired_chips[threshold_index]) {
          ++products[threshold_index].chip_passes[chip];
        }
      }
    }
  }
  std::cerr << "\n";

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;

  // Loop over all thresholds' accumulated products across the full sample.
  for (auto& product : products) {
    std::vector<double> rates_hz;

    // Per threshold, loop through all accumulated chips within.
    for (const auto& [chip, passes] : product.chip_passes) {
      (void)chip;

      // Compute stats.
      const double rate_hz = static_cast<double>(passes) / total_time_sec;
      product.h_rate->Fill(rate_hz);
      rates_hz.push_back(rate_hz);
    }
    if (rates_hz.empty()) continue;
    std::sort(rates_hz.begin(), rates_hz.end());
    product.p95_rate_hz = rates_hz[percentile_index(rates_hz.size(), 0.95)];
    product.p99_rate_hz = rates_hz[percentile_index(rates_hz.size(), 0.99)];
    product.max_rate_hz = rates_hz.back();
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  auto* hist_dir = output.mkdir("hists");
  hist_dir->cd();
  for (const auto& product : products) {
    product.h_rate->Write();
  }

  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    draw_single(output, products[threshold_index], kCoefficients[threshold_index]);
  }
  draw_overlay(output, products);

  output.Close();
  return 0;
}
