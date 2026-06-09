/*

./build/eta_all -i data/reco_bkg_apr -o plots/hit-origins/eta_all.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/MCParticleCollection.h>

#include "classify_hit.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kMCParticleCollection = "MCParticles";

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

double eta(const edm4hep::Vector3d& momentum) {
  const double pt = std::hypot(momentum.x, momentum.y);
  if (pt <= 0.0) return 0.0;
  return std::asinh(momentum.z / pt);
}

void style(TH1D* hist, int color) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetLineWidth(2);
  hist->SetStats(false);
}

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
  for (auto* hist : hists) max = std::max(max, hist->GetMaximum());

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

  const auto args = parse_args(argc, argv);
  const auto files = br::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  std::vector<TH1D*> h_energy_by_origin;
  std::vector<TH1D*> h_eta_by_origin;
  h_energy_by_origin.reserve(kOrigins.size());
  h_eta_by_origin.reserve(kOrigins.size());

  const auto energy_edges = br::log_edges(260, 1e-6, 1e3);
  for (const auto& origin : kOrigins) {
    auto* h_energy = new TH1D(
        (std::string("h_mc_energy_") + origin.label).c_str(),
        "MCParticle origin;E [GeV];Particles",
        260,
        energy_edges.data());
    style(h_energy, origin.color);
    h_energy_by_origin.push_back(h_energy);

    auto* h_eta = new TH1D(
        (std::string("h_mc_eta_") + origin.label).c_str(),
        "MCParticle origin;#eta;Particles",
        180,
        -9.0,
        9.0);
    style(h_eta, origin.color);
    h_eta_by_origin.push_back(h_eta);
  }

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
      if (!br::has_collection(frame, kMCParticleCollection)) continue;
      ++n_events;

      const auto& particles = frame.get<edm4hep::MCParticleCollection>(kMCParticleCollection);
      for (const auto& particle : particles) {
        const int origin = br::origin_index(particle.getGeneratorStatus());
        const double energy = particle.getEnergy();
        h_energy_by_origin[origin]->Fill(energy);
        h_eta_by_origin[origin]->Fill(eta(particle.getMomentum()));
      }
    }
  }

  if (n_events == 0) {
    std::cerr << "No events with " << kMCParticleCollection << " found\n";
    return 1;
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  draw_overlay(output, "MCParticle_energy", "c_mc_energy", "MCParticle origin;E [GeV];Particles", h_energy_by_origin, true, true);
  draw_overlay(output, "MCParticle_eta", "c_mc_eta", "MCParticle origin;#eta;Particles", h_eta_by_origin, false, true);

  std::cout << "\n";
  return 0;
}
