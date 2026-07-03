/*

./build/data_rate_scan -i data/bkg_apr -o plots/data_rates/data_rate_scan.root

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
#include <cstdio>
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
constexpr const char* kHitCollection = "LFHCALHits";
constexpr double kEventWindowSec = 2e-6;
constexpr double kHistMinimum = 0.8;
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;
constexpr std::array<double, 3> kCoefficients = {0.1, 0.5, 1.5};
const std::array<int, kCoefficients.size()> kColors = {kBlue + 1, kBlack, kRed + 1};

struct ThresholdProducts {
  TH1D* h_rate = nullptr;
  // chip_bits[chip]: summed payload bits accumulated across the full sample.
  std::unordered_map<rates::LFHCALChipID, double, rates::LFHCALChipIDHash> chip_bits;
  double p95_rate_gbps = 0.0;
  double p99_rate_gbps = 0.0;
};

struct EventChip {
  // channel_energy[channel]: summed channel signal inside one chip for one event.
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

// Turn a percentile request into the corresponding sorted-data index.
std::size_t percentile_index(std::size_t n_values, double percentile) {
  if (n_values == 0) return 0;
  return static_cast<std::size_t>(percentile * static_cast<double>(n_values - 1));
}


std::string legend_label(double threshold_mip, double percentile_rate_gbps, int percentile) {
  char buffer[160];
  std::snprintf(buffer,
                sizeof(buffer),
                "%s (p%d %.3g Gb/s)",
                threshold_label(threshold_mip).c_str(),
                percentile,
                percentile_rate_gbps);
  return buffer;
}

void draw_single(TDirectory* dir, TH1D* hist, double threshold_mip) {
  dir->cd();

  TCanvas canvas(("c_data_rate_" + threshold_tag(threshold_mip)).c_str(),
                 ("Data rate distribution, threshold " + threshold_label(threshold_mip) + ";data rate [Gb/s];chips").c_str(),
                 1000,
                 800);
  canvas.SetLogx();
  canvas.SetLogy();

  hist->SetStats(false);
  hist->SetLineWidth(2);
  hist->SetMinimum(kHistMinimum);
  hist->Draw("hist");
  canvas.Write();
}

void draw_overlay(TDirectory* dir,
                  const std::array<ThresholdProducts, kCoefficients.size()>& products,
                  int percentile) {
  dir->cd();

  TCanvas canvas(("c_data_rate_threshold_overlay_p" + std::to_string(percentile)).c_str(),
                 "Data rate distributions;data rate [Gb/s];chips",
                 1000,
                 800);
  canvas.SetLogx();
  canvas.SetLogy();

  TLegend legend(0.62, 0.68, 0.88, 0.88);
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
    hist->SetTitle("Data rate distributions;data rate [Gb/s];chips");
    hist->Draw(drew ? "hist same" : "hist");
    drew = true;

    const double percentile_rate_gbps = percentile == 95 ? products[threshold_index].p95_rate_gbps
                                                         : products[threshold_index].p99_rate_gbps;
    legend.AddEntry(hist,
                    legend_label(kCoefficients[threshold_index], percentile_rate_gbps, percentile).c_str(),
                    "l");
  }

  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    const auto& product = products[threshold_index];
    const double percentile_rate_gbps = percentile == 95 ? product.p95_rate_gbps : product.p99_rate_gbps;
    if (percentile_rate_gbps <= 0.0) continue;

    auto* line = new TLine(percentile_rate_gbps, kHistMinimum, percentile_rate_gbps, 1.25 * max_y);
    line->SetLineColor(kColors[threshold_index]);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
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

  const auto rate_edges = rates::log_edges(240, 1.0e-4, 10.0);
  std::array<ThresholdProducts, kCoefficients.size()> products{};
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    const std::string name = "h_data_rate_" + threshold_tag(kCoefficients[threshold_index]);
    const std::string title = "LFHCAL per-chip data rate, threshold " + threshold_label(kCoefficients[threshold_index]) +
        ";data rate [Gb/s];chips";
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

      // Loop all hits in this event.
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        // Increment the corresponding channel's energy for the corresponding readout chip.
        const auto channel = decoder.channel(cell_id);
        auto& event_chip = event_chips[decoder.decode_chip(cell_id)];
        event_chip.channel_energy[channel] += hit.getEnergy();
      }

      // After processing all hits into corresponding channels and chips, loop over all chips.
      for (const auto& [chip, event_chip] : event_chips) {
        std::array<int, kCoefficients.size()> active_counts{};

        // Per chip, loop over all channels.
        for (const auto& [channel, energy_gev] : event_chip.channel_energy) {
          const double mip_gev = rates::mip_energy_gev(channel.rlayerz);

          // Per channel, loop over all thresholds.
          for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {

            // If channel energy meets threshold, count the corresponding payload word.
            if (energy_gev <= kCoefficients[threshold_index] * mip_gev) continue;
            ++active_counts[threshold_index];
          }
        }

        // After summing all fired channels in this chip for this event, convert it into data rate.
        for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
          // Data size = (128 overhead bit + 32 bits per fired channel) * 4 samples
          if ((active_counts[threshold_index]) <= 0) continue;
          const double event_bits =
              (kOverheadBits + kBitsPerHit * static_cast<double>(active_counts[threshold_index])) * kSamplesPerEvent;
          products[threshold_index].chip_bits[chip] += event_bits;
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
    std::vector<double> rates_gbps;

    // Per threshold, loop through all accumulated chips within.
    for (const auto& [chip, total_bits] : product.chip_bits) {
      (void)chip;

      // Compute stats.
      const double rate_gbps = (total_bits / total_time_sec) / 1.0e9;
      product.h_rate->Fill(rate_gbps);
      rates_gbps.push_back(rate_gbps);
    }
    if (rates_gbps.empty()) continue;

    // Sort so percentile lookup becomes a simple indexed read.
    std::sort(rates_gbps.begin(), rates_gbps.end());
    product.p95_rate_gbps = rates_gbps[percentile_index(rates_gbps.size(), 0.95)];
    product.p99_rate_gbps = rates_gbps[percentile_index(rates_gbps.size(), 0.99)];
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

  auto* canvas_dir = output.mkdir("plots");
  for (std::size_t threshold_index = 0; threshold_index < kCoefficients.size(); ++threshold_index) {
    draw_single(canvas_dir, products[threshold_index].h_rate, kCoefficients[threshold_index]);
  }
  draw_overlay(canvas_dir, products, 95);
  draw_overlay(canvas_dir, products, 99);

  output.Close();
  return 0;
}
