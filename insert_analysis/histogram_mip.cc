/*

./build/insert_histogram_mip -i data/mu-_10GeV_insert_100k.edm4hep.root -o insert_plots/histogram_mip.root

4e-4 GeV = MIP
*/

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "utils.h"

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>

namespace fs = std::filesystem;

namespace {

// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------
constexpr const char* kCollectionMatch = "HcalEndcapPInsert";
constexpr int kNBins = 300;
constexpr double kEnergyMinGeV = 0.0;
constexpr double kEnergyMaxGeV = 0.03;

struct Args {
  std::string input_dir;
  std::string output_file;
};

// ----------------------------------------------------------------------------------
// CLI
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
// Helpers
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

void draw_and_write(TDirectory* dir, TH1D* hist, const std::string& canvas_name) {
  // Keep the canvas beside the histogram in the requested output directory.
  dir->cd();
  TCanvas canvas(canvas_name.c_str(), hist->GetTitle(), 1000, 800);
  canvas.SetLogy();
  hist->Draw("hist");
  canvas.Write();
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and expand the input into a file list.
  const auto args = parse_args(argc, argv);
  std::vector<fs::path> files;
  const fs::path input_path(args.input_dir);
  if (fs::is_regular_file(input_path) && input_path.extension() == ".root") {
    files.push_back(input_path);
  } else {
    files = rates::find_root_files(args.input_dir);
  }
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  // Ensure the output directory exists before any event work starts.
  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Create the one histogram that will collect every accepted insert hit.
  auto* hist = new TH1D(
      "h_hit_energy",
      "Insert MIP response;hit energy [GeV];hits",
      kNBins,
      kEnergyMinGeV,
      kEnergyMaxGeV);

  // The decoder is used here as a quick layer filter before filling the spectrum.
  const rates::HcalEndcapPInsertCellIDDecoder decoder;
  std::string collection_name;
  std::uint64_t n_hits = 0;
  rates::FileProgress progress(files.size(), std::cerr);

  // Loop over all files in the input.
  for (const auto& path : files) {
    progress.tick();

    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Per file, loop over all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));

      // Lock onto the first insert hit collection name and reuse it for later events.
      if (collection_name.empty()) {
        collection_name = find_insert_collection(frame);
      }
      if (collection_name.empty() || !rates::has_collection(frame, collection_name)) continue;

      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(collection_name);

      // Loop over all hits in this insert collection.
      for (const auto& hit : hits) {
        const std::uint64_t cell_id = static_cast<std::uint64_t>(hit.getCellID());
        const auto cell = decoder.cell(cell_id);
        if (cell.layer < 1 || cell.layer > 60) continue;
        hist->Fill(hit.getEnergy()); // Add this accepted insert hit to the global energy spectrum.
        ++n_hits; // Track whether the final histogram received any insert hits.
      }
    }
  }
  std::cerr << "\n";

  if (collection_name.empty()) {
    std::cerr << "Failed to find an insert hit collection\n";
    return 1;
  }
  if (n_hits == 0) {
    std::cerr << "No insert hits found in " << args.input_dir << "\n";
    return 1;
  }

  // Write the histogram and its canvas once the full sample has been processed.
  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  hist->Write();
  draw_and_write(&output, hist, "c_hit_energy");

  output.Close();

  return 0;
}
