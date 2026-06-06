/*

./build/energy_res -i data/reco -o plots/performance/energy_res.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ObjectID.h>
#include <podio/ROOTReader.h>

#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>

#include "utils.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>

namespace fs = std::filesystem;

namespace {

constexpr const char* kClusterCollection = "LFHCALClusters";
constexpr const char* kClusterAssocCollection = "LFHCALClusterAssociations";
constexpr double kDefaultMinRecoGeV = 1e-3;
constexpr double kDefaultMinPurity = 0;
constexpr double kDefaultTruthMinGeV = 0.1;       // Used to set the min x-axis range.
constexpr double kDefaultTruthMaxGeV = 1000.0;    // Max x-axis range.
constexpr double kDefaultResponseAbsMax = 3.0;    // Sets the y range.

struct Args {
  std::string input_dir;
  std::string output_file;
  double min_reco_geV = kDefaultMinRecoGeV;
  double min_purity = kDefaultMinPurity;
  double truth_min_geV = kDefaultTruthMinGeV;
  double truth_max_geV = kDefaultTruthMaxGeV;
  double response_abs_max = kDefaultResponseAbsMax;
};

// ----------------------------- handle CLI inputs -----------------------------
void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0
            << " -i INPUT_DIR -o OUTPUT.root [-m MIN_RECO_GEV] [-p MIN_PURITY]"
            << " [--truth-min GEV] [--truth-max GEV] [--response-max VAL]\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_dir = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      args.output_file = argv[++i];
    } else if ((arg == "-m" || arg == "--min-reco") && i + 1 < argc) {
      args.min_reco_geV = std::stod(argv[++i]);
    } else if ((arg == "-p" || arg == "--min-purity") && i + 1 < argc) {
      args.min_purity = std::stod(argv[++i]);
    } else if (arg == "--truth-min" && i + 1 < argc) {
      args.truth_min_geV = std::stod(argv[++i]);
    } else if (arg == "--truth-max" && i + 1 < argc) {
      args.truth_max_geV = std::stod(argv[++i]);
    } else if (arg == "--response-max" && i + 1 < argc) {
      args.response_abs_max = std::stod(argv[++i]);
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

std::uint64_t object_key(const podio::ObjectID& id) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(id.collectionID)) << 32) |
         static_cast<std::uint32_t>(id.index);
}

void draw_and_write(TDirectory* dir, TH1* hist, const char* canvas_name, bool logz = false, bool logx = false) {
  dir->cd();
  TCanvas canvas(canvas_name, hist->GetTitle(), 1000, 800);
  if (logz) canvas.SetLogz();
  if (logx) canvas.SetLogx();
  hist->SetStats(false);
  hist->Draw(hist->InheritsFrom(TH2::Class()) ? "colz" : "hist");
  canvas.Write();
}

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and find reco event files.
  const auto args = parse_args(argc, argv);
  const auto files = br::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  // Ensure output.
  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Compute bin edges that's evenly spaced in log scale.
  const auto truth_edges = br::log_edges(200, args.truth_min_geV, args.truth_max_geV);

  // Initialize TH2D.
  auto* h_truth_vs_response = new TH2D(
      "h_truth_vs_response",
      "LFHCAL cluster response;E_{truth}^{matched} [GeV];(E_{reco} - E_{truth}^{matched}) / E_{truth}^{matched}",
      200, truth_edges.data(), // .data() gives a pointer starting at the first element of the vector’s internal array
      200, -args.response_abs_max, args.response_abs_max);

  std::uint64_t n_events = 0;
  std::uint64_t n_clusters_filled = 0;

  br::FileProgress progress(files.size(), std::cerr);

  // Loop through all files.
  for (const auto& path : files) {
    progress.tick();

    // Grab the events for this file.
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Loop through all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));

      // Skip this event if collections don't exist.
      if (!br::has_collection(frame, kClusterCollection) || !br::has_collection(frame, kClusterAssocCollection)) {
        continue;
      }

      // Otherwise, grab collections.
      ++n_events;
      const auto& clusters = frame.get<edm4eic::ClusterCollection>(kClusterCollection);
      const auto& assocs = frame.get<edm4eic::MCRecoClusterParticleAssociationCollection>(kClusterAssocCollection);

      std::unordered_map<std::uint64_t, edm4eic::MCRecoClusterParticleAssociation> best_assoc_by_cluster;
      best_assoc_by_cluster.reserve(assocs.size());

      // Loop through all associations.
      for (const auto& assoc : assocs) {

        // For each associatioon, build a (key, value) map:

        //  key = reco cluster ID
        const auto key = object_key(assoc.getRec().getObjectID());
        auto it = best_assoc_by_cluster.find(key);

        //  value = highest-weight association seen for that cluster.
        if (it == best_assoc_by_cluster.end() || assoc.getWeight() > it->second.getWeight()) {
          best_assoc_by_cluster[key] = assoc;
        }
      }

      // Loop through all clusters.
      for (const auto& cluster : clusters) {
        const double reco_energy = cluster.getEnergy();
        if (reco_energy <= args.min_reco_geV) continue; // Filter by min energy.

        // Grab this cluster's best association from the map.
        const auto assoc_it = best_assoc_by_cluster.find(object_key(cluster.getObjectID()));
        if (assoc_it == best_assoc_by_cluster.end()) continue;  // Skip clusters w/ no truth associations.
        const auto& best_assoc = assoc_it->second;
        if (best_assoc.getWeight() < args.min_purity) continue; // Skip clusters with purity below the user input.

        // Require valid truth energy.
        const double truth_energy = best_assoc.getSim().getEnergy();
        if (truth_energy <= 0.0) continue;

        // Fill histogram.
        h_truth_vs_response->Fill(truth_energy, (reco_energy - truth_energy) / truth_energy);
        ++n_clusters_filled;
      }
    }
  }

  if (n_events == 0) {
    std::cerr << "No events with " << kClusterCollection << " and " << kClusterAssocCollection
              << " found\n";
    return 1;
  }
  if (n_clusters_filled == 0) {
    std::cerr << "No truth-matched clusters found\n";
    return 1;
  }

  // Write outputs.
  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  auto* hist_dir = output.mkdir("hists");
  hist_dir->cd();
  h_truth_vs_response->Write();

  auto* canvas_dir = output.mkdir("canvases");
  draw_and_write(canvas_dir, h_truth_vs_response, "c_truth_vs_response", true, true);

  std::cout << "\n";
  return 0;
}
