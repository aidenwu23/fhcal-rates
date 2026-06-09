/*

./build/hit_origins -i data/reco_bkg_feb -o plots/hit-origins/origins.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TLegend.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "classify_hit.h"
#include "decode_channel.h"
#include "edep.h"
#include "eta.h"
#include "utils.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

namespace {

const std::string kHitCollection = "LFHCALHits";
constexpr int kNReadoutLayers = 7;
constexpr double kEventWindowSec = 2e-6;
const std::vector<double> kEtaThresholdsGeV = {0.0, 5e-4};  // Min energy threshold applied to hits in the eta graphs.

struct OriginInfo {
  const char* label;
  int color;
};

const std::vector<OriginInfo> kOrigins = {
    {"signal", kBlack},
    {"synrad", kRed + 1},
    {"eBrem", kBlue + 1},
    {"eTouschek", kMagenta + 1},
    {"eCoulomb", kGreen + 2},
    {"pBeamGas", kOrange + 7},
};

// ----------------------------- handle CLI inputs -----------------------------
struct Args {
  std::string input_dir;
  std::string output_file;
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
// -----------------------------------------------------------------------------

// Draw an overlay canvas by looping over all provided hists.
void draw_overlay(TFile& file,
                  const char* directory,
                  const char* canvas_name,
                  const char* title,
                  const std::vector<TH1D*>& hists,
                  bool logx,
                  bool logy) {
  file.mkdir(directory)->cd();

  TCanvas canvas(canvas_name, title, 1000, 800);
  if (logx) canvas.SetLogx();
  if (logy) canvas.SetLogy();

  TLegend legend(0.66, 0.62, 0.88, 0.88);
  legend.SetBorderSize(1);

  double max = 0.0;

  for (auto* hist : hists) {
    max = std::max(max, hist->GetMaximum());
  }

  for (std::size_t i = 0; i < hists.size(); ++i) {
    auto* hist = hists[i];
    hist->Write();
    hist->SetTitle(title);
    if (logy) {
      hist->SetMinimum(0.5);
      if (max > 0.0) hist->SetMaximum(max * 5.0);
    }
    hist->Draw(i == 0 ? "hist" : "hist same");
    legend.AddEntry(hist, kOrigins[i].label, "l");
  }

  legend.Draw();
  canvas.Write();
  file.cd();
}

// Write a per-origin hit-rate summary as CSV.
void write_rate_csv(const fs::path& output_path,
                    const std::vector<std::uint64_t>& origin_hit_counts,
                    std::uint64_t n_events) {
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::ofstream out(output_path);
  if (!out) {
    throw std::runtime_error("Failed to open output CSV");
  }

  out << "origin_label,total_hits,rate_hz\n";
  out << std::setprecision(10);

  for (std::size_t i = 0; i < kOrigins.size(); ++i) {
    const double rate_hz = total_time_sec > 0.0
                               ? static_cast<double>(origin_hit_counts[i]) / total_time_sec // Hz
                               : 0.0;
    out << kOrigins[i].label << ','
        << origin_hit_counts[i] << ','
        << rate_hz << '\n';
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args.
  const auto args = parse_args(argc, argv);
  const auto& input_dir = args.input_dir;
  const auto& output_file = args.output_file;

  // Find input files.
  const auto files = br::find_root_files(input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << input_dir << "\n";
    return 1;
  }

  // Handle output.
  fs::path output_path = output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  TFile output(output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << output_file << "\n";
    return 1;
  }

  // Make a decoder.
  const br::LFHCALDecoder decoder;
  const auto edges = br::log_edges(260, 1e-10, 10.0); // Log edges for nicer display.

  // Used for styling histograms.
  std::vector<const char*> origin_labels;
  std::vector<int> origin_colors;

  origin_labels.reserve(kOrigins.size());
  origin_colors.reserve(kOrigins.size());

  // Extract the label and color fields from kOrigins into two parallel vectors so the helper functions can use them.
  for (const auto& origin : kOrigins) {
    origin_labels.push_back(origin.label);
    origin_colors.push_back(origin.color);
  }

  // Make a histogram for all 7 readout layers.
  std::vector<br::origins::LayerEdepHists> edep_layers;
  edep_layers.reserve(kNReadoutLayers);
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    // For each layer, make a histogram for all generatorStatus families.
    edep_layers.push_back(br::origins::make_edep_hists(origin_labels, origin_colors, layer, 260, edges.data()));
  }

  // Also make an inclusive version for all 7.
  auto edep_hit = br::origins::make_summed_edep_hists(origin_labels, origin_colors, 260, edges.data());

  // Make a histogram for eta hists.
  std::vector<br::origins::ThresholdEtaHists> eta_hists;
  eta_hists.reserve(kEtaThresholdsGeV.size());

  // One of the thresholds is 0 so its just <no threshold, threshold>.
  for (double threshold_geV : kEtaThresholdsGeV) { 
    eta_hists.push_back(br::origins::make_eta_hists(origin_labels, origin_colors, threshold_geV));
  }

  TH1I status("h_status", "Raw generatorStatus;generatorStatus;Contributions", 8000, -1000, 7000);
  status.SetStats(false);
  std::vector<std::uint64_t> origin_hit_counts(kOrigins.size(), 0);
  std::uint64_t n_events = 0;

  br::FileProgress progress(files.size(), std::cerr);

  // Loop through files.
  for (const auto& path : files) {
    progress.tick();

    // Read file.
    podio::ROOTReader reader;
    reader.openFile(path.string());

    // Loop through all events.
    const std::size_t total_events = reader.getEntries("events");
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {

      auto data = reader.readEntry("events", event_index);
      if (!data) {
        continue;
      }

      podio::Frame frame(std::move(data));
      if (!br::has_collection(frame, kHitCollection)) {
        continue;
      }
      ++n_events;

      // Loop through all hits for this event.
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);
      for (const auto& hit : hits) {
        const int layer = decoder.get(static_cast<std::uint64_t>(hit.getCellID()), "rlayerz");

        // Per hit, loop through all contributions.
        for (const auto& contribution : hit.getContributions()) {
          const auto particle = contribution.getParticle();
          const int generator_status = particle.getGeneratorStatus();

          // Convert generatorStatus one of 6 origin indices.
          // 0 = signal
          // 1 = synrad
          // 2 = eBrem
          // 3 = eTouschek
          // 4 = eCoulomb
          // 5 = pBeamGas
          const int contribution_origin = br::origin_index(generator_status);
          const double contribution_energy = contribution.getEnergy();

          // Record raw generator status of the contribution (for debugging)
          status.Fill(generator_status);
          if (contribution_energy <= 0.0) continue;

          ++origin_hit_counts[contribution_origin];

          if (layer >= 0 && layer < kNReadoutLayers) {
            br::origins::fill_edep(edep_layers[layer], contribution_origin, contribution_energy);
          }

          for (auto& eta_hists_for_threshold : eta_hists) {
            br::origins::fill_eta(eta_hists_for_threshold, contribution_origin, hit.getPosition(), contribution_energy);
          }
        }
      }
    }
  }

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  for (const auto& layer_hists : edep_layers) {
    br::origins::sum_edep_into(edep_hit, layer_hists);
  }

  draw_overlay(output, "Edep_hit", "c_edep_hit", "LFHCAL hit origin;E_{dep} [GeV];Hits", edep_hit, true, true);
  for (const auto& layer_hists : edep_layers) {
    const auto layer_dir = "Edep_hit_layer" + std::to_string(layer_hists.layer);
    const auto layer_canvas = "c_edep_hit_layer" + std::to_string(layer_hists.layer);
    draw_overlay(output,
                 layer_dir.c_str(),
                 layer_canvas.c_str(),
                 layer_hists.title.c_str(),
                 layer_hists.by_origin,
                 true,
                 true);
  }
  for (const auto& eta_hists_for_threshold : eta_hists) {
    draw_overlay(output,
                 eta_hists_for_threshold.directory.c_str(),
                 eta_hists_for_threshold.canvas_name.c_str(),
                 eta_hists_for_threshold.title.c_str(),
                 eta_hists_for_threshold.by_origin,
                 false,
                 true);
  }
  output.mkdir("Status")->cd();
  status.Write();
  output.Close();

  try {
    auto csv_path = output_path;
    csv_path.replace_extension(".csv");
    write_rate_csv(csv_path, origin_hit_counts, n_events);
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << '\n';
    return 1;
  }

  std::cout << "\n";
  return 0;
}
