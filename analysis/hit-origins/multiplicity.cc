/*

./build/multiplicity -i data/bkg_apr -o plots/hit-origins/multiplicity.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1I.h>

#include <podio/Frame.h>
#include <podio/ObjectID.h>
#include <podio/ROOTReader.h>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "classify_hit.h"
#include "utils.h"

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_set>

namespace fs = std::filesystem;

namespace {

constexpr const char* kHitCollection = "LFHCALHits";

// ----------------------------------------------------------------------------------
// CLI handling.
// ----------------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------------
// Helpers.
// ----------------------------------------------------------------------------------
std::uint64_t object_key(const podio::ObjectID& id) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(id.collectionID)) << 32) |
         static_cast<std::uint32_t>(id.index);
}

void draw_and_write(TDirectory* dir, TH1* hist, const char* canvas_name) {
  dir->cd();
  TCanvas canvas(canvas_name, hist->GetTitle(), 1000, 800);
  hist->SetStats(false);
  hist->Draw("hist");
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

  TH1I h_n_unique_contributors("h_n_unique_contributors",
                               "LFHCAL hit contributor multiplicity;Unique contributing particles per hit;Hits",
                               100, 0, 100);
  TH1I h_n_unique_families("h_n_unique_families",
                           "LFHCAL hit origin-family multiplicity;Unique origin families per hit;Hits",
                           10, 0, 10);

  std::uint64_t n_events = 0;
  br::FileProgress progress(files.size(), std::cerr);

  // Loop files.
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

      // Loop hits.
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);
      for (const auto& hit : hits) {
        std::unordered_set<std::uint64_t> unique_contributors;
        std::unordered_set<int> unique_families;

        // Loop contributions and insert their ObjectID + generatorStatus family into unordered sets.
        for (const auto& contribution : hit.getContributions()) {
          unique_contributors.insert(object_key(contribution.getParticle().getObjectID()));
          unique_families.insert(br::origin_index(contribution.getParticle().getGeneratorStatus()));
        }

        // Size of each unordered set = number of unique elements.
        h_n_unique_contributors.Fill(static_cast<int>(unique_contributors.size()));
        h_n_unique_families.Fill(static_cast<int>(unique_families.size()));
      }
    }
  }

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  auto* hist_dir = output.mkdir("hists");
  hist_dir->cd();
  h_n_unique_contributors.Write();
  h_n_unique_families.Write();

  auto* canvas_dir = output.mkdir("canvases");
  draw_and_write(canvas_dir, &h_n_unique_contributors, "c_n_unique_contributors");
  draw_and_write(canvas_dir, &h_n_unique_families, "c_n_unique_families");

  std::cout << "\n";
  return 0;
}
