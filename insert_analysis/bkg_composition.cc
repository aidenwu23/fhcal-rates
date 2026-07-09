/*

./build/bkg_composition -i data/bkg_apr -o insert_plots/bkg_composition.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "classify_hit.h"
#include "utils.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

namespace {

// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
constexpr const char* kCollectionMatch = "HcalEndcapPInsert";

struct OriginInfo {
  const char* label;
  int color;
};

const std::vector<OriginInfo> kOrigins = {
    {"DIS", kBlack},
    {"synrad", kRed + 1},
    {"eBrem", kBlue + 1},
    {"eTouschek", kMagenta + 1},
    {"eCoulomb", kGreen + 2},
    {"pBeamGas", kOrange + 7},
    {"other", kGray + 2},
};

struct Args {
  std::string input_dir;
  std::string output_file;
};

// ----------------------------------------------------------------------------------
// CLI handling.
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
// Helpers and plotters.
// ----------------------------------------------------------------------------------
std::string find_insert_collection(const podio::Frame& frame) {
  for (const auto& name : frame.getAvailableCollections()) {
    if (name.find(kCollectionMatch) != std::string::npos &&
        name.find("Contributions") == std::string::npos) {
      return name;
    }
  }
  return "";
}

void style(TH1D* hist, int color) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetLineWidth(2);
  hist->SetStats(false);
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main.
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  const auto args = parse_args(argc, argv);
  const auto files = rates::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  std::vector<TH1D*> hists;
  hists.reserve(kOrigins.size());
  const auto edges = rates::log_edges(260, 1e-10, 10.0);
  for (const auto& origin : kOrigins) {
    const std::string name = std::string("h_edep_hit_") + origin.label;
    auto* hist = new TH1D(name.c_str(), "Insert hit origin;E_{dep} [GeV];Hits", 260, edges.data());
    style(hist, origin.color);
    hists.push_back(hist);
  }

  std::string insert_collection_name;
  std::uint64_t n_events = 0;
  rates::FileProgress progress(files.size(), std::cerr);

  // Loop thru all files.
  for (const auto& path : files) {
    progress.tick();

    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Per file, loop thru all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (insert_collection_name.empty()) {
        insert_collection_name = find_insert_collection(frame);
      }
      if (insert_collection_name.empty() || !rates::has_collection(frame, insert_collection_name)) continue;
      ++n_events;

      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(insert_collection_name);

      // Per event, loop thru all hits.
      for (const auto& hit : hits) {
        std::vector<double> energy_by_origin(kOrigins.size(), 0.0);

        // Per hit, loop thru all hit contributioons.
        for (const auto& contribution : hit.getContributions()) {

          // Convert generatorStatus into an origin type and increment the energy.
          const int origin = rates::origin_index(contribution.getParticle().getGeneratorStatus());
          const double energy = contribution.getEnergy();
          if (origin < 0 || energy <= 0.0) continue;
          energy_by_origin[origin] += energy;
        }

        // After accumulating all hit contributions by origin, histogram this hit.
        for (std::size_t origin = 0; origin < energy_by_origin.size(); ++origin) {
          if (energy_by_origin[origin] <= 0.0) continue;
          hists[origin]->Fill(energy_by_origin[origin]);
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

  output.cd();
  TCanvas canvas("c_edep_hit", "Insert hit origin", 1000, 800);
  canvas.SetLogx();
  canvas.SetLogy();

  TLegend legend(0.66, 0.62, 0.88, 0.88);
  legend.SetBorderSize(1);

  double max_value = 0.0;
  for (auto* hist : hists) {
    max_value = std::max(max_value, hist->GetMaximum());
  }

  for (std::size_t i = 0; i < hists.size(); ++i) {
    auto* hist = hists[i];
    hist->SetMinimum(0.5);
    if (max_value > 0.0) hist->SetMaximum(max_value * 5.0);
    hist->Write();
    hist->Draw(i == 0 ? "hist" : "hist same");
    legend.AddEntry(hist, kOrigins[i].label, "l");
  }

  legend.Draw();
  canvas.Write();
  output.Close();
  return 0;
}
