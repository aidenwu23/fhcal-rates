/*

./build/hit_origins -i data/reco -o plots/hit-origins/origins.root

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
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

const std::string kHitCollection = "LFHCALHits";

struct OriginInfo {
  const char* key;
  const char* label;
  int color;
};

// Color mapping.
const std::vector<OriginInfo> kOrigins = {
    {"signal", "Signal", kBlack},
    {"synrad", "Synrad", kRed + 1},
    {"ebrem", "eBrem", kBlue + 1},
    {"etouschek", "eTouschek", kMagenta + 1},
    {"ecoulomb", "eCoulomb", kGreen + 2},
    {"pbeamgas", "pBeamGas", kOrange + 7},
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

// Compute eta.
double eta(const edm4hep::Vector3f& position) {
  const double pt = std::hypot(position.x, position.y);
  if (pt <= 0.0) return 0.0;
  return std::asinh(position.z / pt);
}

// Build log spaced bin edges (so things display nicely on log scale).
std::vector<double> log_edges(int bins, double low, double high) {
  std::vector<double> edges(bins + 1);
  const double log_low = std::log10(low);
  const double log_high = std::log10(high);
  for (int bin = 0; bin <= bins; ++bin) {
    const double fraction = static_cast<double>(bin) / bins;
    edges[bin] = std::pow(10.0, log_low + fraction * (log_high - log_low));
  }
  return edges;
}

// Style a histogram.
void style(TH1D* hist, int color) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetLineWidth(2);
  hist->SetStats(false);
}

// Draw an overlay canvas by looping over all hists.
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

  const auto edges = log_edges(260, 1e-10, 10.0);

  // Create two histogram types. One bins by energy deposited (edep) and the other bins by eta.
  std::vector<TH1D*> edep_hit;
  std::vector<TH1D*> eta_hit;
  for (const auto& origin : kOrigins) {
    auto* hit_energy = new TH1D((std::string("h_edep_hit_") + origin.key).c_str(),
                                "LFHCAL hit origin;E_{dep} [GeV];Hits", 260, edges.data());
    auto* hit_eta = new TH1D((std::string("h_eta_hit_") + origin.key).c_str(),
                             "LFHCAL hit origin;#eta;Hits", 240, -12.0, 12.0);
    style(hit_energy, origin.color);
    style(hit_eta, origin.color);
    edep_hit.push_back(hit_energy);
    eta_hit.push_back(hit_eta);
  }

  TH1I status("h_status", "Raw generatorStatus;generatorStatus;Contributions", 8000, -1000, 7000);
  status.SetStats(false);

  const auto total_files = files.size();
  std::size_t i = 0;

  // Loop through files.
  for (const auto& path : files) {
    ++i;
    std::cerr << "\r" << i << "/" << total_files << " files read." << std::flush;

    // Read file.
    podio::ROOTReader reader;
    reader.openFile(path.string());

    // Loop through all events.
    const std::size_t total_events = reader.getEntries("events");
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {

      // Grab collection.
      auto data = reader.readEntry("events", event_index);
      if (!data) {
        continue;
      }

      podio::Frame frame(std::move(data));
      if (!br::has_collection(frame, kHitCollection)) {
        continue;
      }

      // Loop through all hits for this event.
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);
      for (const auto& hit : hits) {
        double biggest_contribution = -1.0;
        int hit_origin = 0;

        // Per hit loop through all contributions.
        for (const auto& contribution : hit.getContributions()) {
          const auto particle = contribution.getParticle();
          const int generator_status = particle.getGeneratorStatus();
          const int contribution_origin = br::origin_index(generator_status);
          const double contribution_energy = contribution.getEnergy();

          // Record generator status of the contribution.
          status.Fill(generator_status);

          // Define a hit origin by the dominant energy contributor.
          if (contribution_energy > biggest_contribution) {
            biggest_contribution = contribution_energy;
            hit_origin = contribution_origin;
          }
        }

        // Fill hists.
        if (hit.getEnergy() > 0.0) {
          edep_hit[hit_origin]->Fill(hit.getEnergy());
        }
        eta_hit[hit_origin]->Fill(eta(hit.getPosition()));
      }
    }
  }

  draw_overlay(output, "Edep_hit", "c_edep_hit", "LFHCAL hit origin;E_{dep} [GeV];Hits", edep_hit, true, true);
  draw_overlay(output, "Eta_hit", "c_eta_hit", "LFHCAL hit origin;#eta;Hits", eta_hit, false, true);
  output.mkdir("Status")->cd();
  status.Write();
  output.Close();

  std::cout << "\n";
  return 0;
}
