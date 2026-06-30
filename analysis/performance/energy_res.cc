/*

./build/energy_res -i data/bkg_apr -o plots/performance/energy_res.root
./build/energy_res -i data/reco_no_bkg_apr -o plots/performance/energy_res.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ObjectID.h>
#include <podio/ROOTReader.h>

#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>

#include "utils.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kClusterCollection = "LFHCALClusters";
constexpr const char* kClusterAssocCollection = "LFHCALClusterAssociations";
constexpr double kDefaultMinRecoGeV = 1e-3;
constexpr double kDefaultMinPurity = 0;
constexpr double kDefaultTruthMinGeV = 0.1;       // Used to set the min x-axis range.
constexpr double kDefaultTruthMaxGeV = 1000.0;    // Max x-axis range.
constexpr double kDefaultResponseAbsMax = 3.0;    // Sets the y range.

struct RatioRange {
  double low_geV;
  double high_geV;
  const char* name;
  const char* title;
};

constexpr RatioRange kRatioRanges[] = {
    {0.0, 0.0, "h_ratio_inclusive", "LFHCAL cluster energy ratio, inclusive;E_{reco}^{cluster} / E_{truth}^{matched};Clusters"},
    {0.0, 10.0, "h_ratio_0_10", "LFHCAL cluster energy ratio;E_{reco}^{cluster} / E_{truth}^{matched};Clusters"},
    {10.0, 20.0, "h_ratio_10_20", "LFHCAL cluster energy ratio;E_{reco}^{cluster} / E_{truth}^{matched};Clusters"},
    {20.0, 40.0, "h_ratio_20_40", "LFHCAL cluster energy ratio;E_{reco}^{cluster} / E_{truth}^{matched};Clusters"},
    {40.0, 80.0, "h_ratio_40_80", "LFHCAL cluster energy ratio;E_{reco}^{cluster} / E_{truth}^{matched};Clusters"},
};

struct Args {
  std::string input_dir;
  std::string output_file;
  double min_reco_geV = kDefaultMinRecoGeV;
  double min_purity = kDefaultMinPurity;
  double truth_min_geV = kDefaultTruthMinGeV;
  double truth_max_geV = kDefaultTruthMaxGeV;
  double response_abs_max = kDefaultResponseAbsMax;
};

struct Counters {
  std::uint64_t n_events = 0;
  std::uint64_t n_clusters_filled = 0;
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

// For the 2D histogram.
void draw_and_write(TDirectory* dir, TH1* hist, const char* canvas_name, bool logz = false, bool logx = false) {
  dir->cd();
  TCanvas canvas(canvas_name, hist->GetTitle(), 1000, 800);
  if (logz) canvas.SetLogz();
  if (logx) canvas.SetLogx();
  hist->SetStats(false);
  hist->Draw(hist->InheritsFrom(TH2::Class()) ? "colz" : "hist");
  canvas.Write();
}

// For the other histogram type.
void draw_and_write(TDirectory* dir, TGraphErrors* graph, const char* canvas_name, bool logx = false) {
  dir->cd();
  TCanvas canvas(canvas_name, graph->GetTitle(), 1000, 800);
  if (logx) canvas.SetLogx();
  graph->SetMarkerStyle(20);
  graph->SetLineWidth(2);
  graph->Draw("AP");
  canvas.Write();
}

// --------------------------------- Observables ---------------------------------
void fill_truth_vs_response(double reco_energy, double truth_energy, TH2D* h_truth_vs_response) {
    h_truth_vs_response->Fill(truth_energy, (reco_energy - truth_energy) / truth_energy);
}

// Output is fed into finalize_resolution.
void fill_truth_vs_reco(double reco_energy, double truth_energy, TH2D* h_truth_vs_reco) {
  h_truth_vs_reco->Fill(truth_energy, reco_energy);
}

// "finalize" the truth vs reco histogram by adding error bars for each bin once all data has been collected.
void finalize_resolution(TH2D* h_truth_vs_reco, TGraphErrors* g_truth_vs_resolution) {
  int resolution_point = 0;

  for (int x_bin = 1; x_bin <= h_truth_vs_reco->GetNbinsX(); ++x_bin) {
    auto proj_name = std::string("h_reco_bin_") + std::to_string(x_bin);
    
    // Take "h_truth_vs_reco" (hist) and extract the y values (E_reco) for one x bin (E_truth).
    auto proj = std::unique_ptr<TH1D>(h_truth_vs_reco->ProjectionY(proj_name.c_str(), x_bin, x_bin));
    
    if (proj->GetEntries() < 2) continue;

    const double low_edge = h_truth_vs_reco->GetXaxis()->GetBinLowEdge(x_bin);
    const double up_edge = h_truth_vs_reco->GetXaxis()->GetBinUpEdge(x_bin);
    const double truth_center = std::sqrt(low_edge * up_edge);
    const double sigma_reco = proj->GetStdDev();
    const double sigma_reco_err = proj->GetStdDevError();

    g_truth_vs_resolution->SetPoint(resolution_point, truth_center, sigma_reco / truth_center);
    g_truth_vs_resolution->SetPointError(resolution_point, 0.0, sigma_reco_err / truth_center);
    ++resolution_point;
  }
}

// Fill a reco_energy/truth_energy histogram for multiple truth energy GeV ranges.
// 0-10 GeV, 10-20 GeV, 20-40 GeV, 40-80 GeV, and inclusive.
void fill_ratio_histograms(double reco_energy, double truth_energy, const std::vector<TH1D*>& h_ratio_hists) {
  const double energy_ratio = reco_energy / truth_energy;
  h_ratio_hists[0]->Fill(energy_ratio);

  for (std::size_t i = 1; i < std::size(kRatioRanges); ++i) {
    const auto& range = kRatioRanges[i];
    if (truth_energy >= range.low_geV && truth_energy < range.high_geV) {
      h_ratio_hists[i]->Fill(energy_ratio);
      break;
    }
  }
}
// -----------------------------------------------------------------------------

}  // namespace

int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and find reco event files.
  const auto args = parse_args(argc, argv);
  const auto files = rates::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  // Ensure output.
  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Compute bin edges that's evenly spaced in log scale.
  const auto truth_edges = rates::log_edges(200, args.truth_min_geV, args.truth_max_geV);

  // Initialize graphs.
  auto* h_truth_vs_response = new TH2D(
      "h_truth_vs_response",
      "LFHCAL cluster response;E_{truth}^{matched} [GeV];(E_{reco}^{cluster} - E_{truth}^{matched}) / E_{truth}^{matched}",
      200, truth_edges.data(), // .data() gives a pointer starting at the first element of the vector’s internal array
      200, -args.response_abs_max, args.response_abs_max);

  auto* h_truth_vs_reco = new TH2D(
      "h_truth_vs_reco",
      "LFHCAL reco energy;E_{truth}^{matched} [GeV];E_{reco}^{cluster} [GeV]",
      200, truth_edges.data(),
      200, truth_edges.data());

  auto* g_truth_vs_resolution = new TGraphErrors();
  g_truth_vs_resolution->SetName("g_truth_vs_resolution");
  g_truth_vs_resolution->SetTitle("LFHCAL energy resolution;E_{truth}^{matched} [GeV];#sigma(E_{reco}^{cluster}) / E_{truth}^{matched}");

  std::vector<TH1D*> h_ratio_hists;
  h_ratio_hists.reserve(std::size(kRatioRanges));
  for (const auto& range : kRatioRanges) {
    h_ratio_hists.push_back(new TH1D(range.name, range.title, 120, 0.0, 3.0));
  }

  Counters counters;
  rates::FileProgress progress(files.size(), std::cerr);

  // Loop thru all input files.
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
      if (!rates::has_collection(frame, kClusterCollection) || !rates::has_collection(frame, kClusterAssocCollection)) {
        continue;
      }

      ++counters.n_events;

      // Find the best truth association/particle to each cluster.
      const auto& clusters = frame.get<edm4eic::ClusterCollection>(kClusterCollection);
      const auto& assocs = frame.get<edm4eic::MCRecoClusterParticleAssociationCollection>(kClusterAssocCollection);
      std::unordered_map<std::uint64_t, edm4eic::MCRecoClusterParticleAssociation> best_assoc_by_cluster;
      best_assoc_by_cluster.reserve(assocs.size());

      for (const auto& assoc : assocs) {
        const auto key = object_key(assoc.getRec().getObjectID());
        auto it = best_assoc_by_cluster.find(key);
        if (it == best_assoc_by_cluster.end() || assoc.getWeight() > it->second.getWeight()) {
          best_assoc_by_cluster[key] = assoc;
        }
      }

      // For this event, loop through all clusters (E_reco).
      for (const auto& cluster : clusters) {
        const double reco_energy = cluster.getEnergy();
        if (reco_energy <= args.min_reco_geV) continue;

        const auto assoc_it = best_assoc_by_cluster.find(object_key(cluster.getObjectID()));
        if (assoc_it == best_assoc_by_cluster.end()) continue;
        const auto& best_assoc = assoc_it->second;  // Associated truth energy to this cluster.
        if (best_assoc.getWeight() < args.min_purity) continue;

        const double truth_energy = best_assoc.getSim().getEnergy();
        if (truth_energy <= 0.0) continue;

        fill_truth_vs_response(reco_energy, truth_energy, h_truth_vs_response);
        fill_truth_vs_reco(reco_energy, truth_energy, h_truth_vs_reco);
        fill_ratio_histograms(reco_energy, truth_energy, h_ratio_hists);
        ++counters.n_clusters_filled;
      }
    }
  }

  if (counters.n_events == 0) {
    std::cerr << "No events with " << kClusterCollection << " and " << kClusterAssocCollection
              << " found\n";
    return 1;
  }
  if (counters.n_clusters_filled == 0) {
    std::cerr << "No truth-matched clusters found\n";
    return 1;
  }

  finalize_resolution(h_truth_vs_reco, g_truth_vs_resolution);

  // Write outputs.
  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  auto* hist_dir = output.mkdir("hists");
  hist_dir->cd();
  h_truth_vs_response->Write();
  h_truth_vs_reco->Write();
  g_truth_vs_resolution->Write();
  for (auto* hist : h_ratio_hists) hist->Write();

  auto* canvas_dir = output.mkdir("canvases");
  draw_and_write(canvas_dir, h_truth_vs_response, "c_truth_vs_response", true, true);
  draw_and_write(canvas_dir, g_truth_vs_resolution, "c_truth_vs_resolution", true);
  draw_and_write(canvas_dir, h_ratio_hists[0], "c_ratio_inclusive");
  draw_and_write(canvas_dir, h_ratio_hists[1], "c_ratio_0_10");
  draw_and_write(canvas_dir, h_ratio_hists[2], "c_ratio_10_20");
  draw_and_write(canvas_dir, h_ratio_hists[3], "c_ratio_20_40");
  draw_and_write(canvas_dir, h_ratio_hists[4], "c_ratio_40_80");

  std::cout << "\n";
  return 0;
}
