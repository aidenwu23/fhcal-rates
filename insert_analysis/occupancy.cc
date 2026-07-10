/*

./build/insert_occupancy -i data/bkg_apr -o insert_plots/occupancy.root

*/

#include <TFile.h>
#include <TH1.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "event_hit.h"
#include "insert_to_lfhcal.h"
#include "lfhcal_segments/include/no_inner_ring.h"
#include "lfhcal_segments/include/with_inner_ring.h"
#include "physical_layers/include/no_inner_ring.h"
#include "physical_layers/include/with_inner_ring.h"
#include "utils.h"

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
constexpr const char* kInsertCollectionMatch = "HcalEndcapPInsert";

struct Args {
  std::string input_path;
  std::string output_file;
};

// ----------------------------------------------------------------------------------
// CLI handling.
// ----------------------------------------------------------------------------------
void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT.root|DIR -o OUTPUT.root\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_path = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      args.output_file = argv[++i];
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input_path.empty() || args.output_file.empty()) {
    usage(argv[0]);
    std::exit(1);
  }

  return args;
}

// ----------------------------------------------------------------------------------
// Helper(s).
// ----------------------------------------------------------------------------------
std::string find_insert_collection(const podio::Frame& frame) {
  // Contributions are stored separately from the simulated hit collection.
  for (const auto& name : frame.getAvailableCollections()) {
    if (name.find(kInsertCollectionMatch) != std::string::npos &&
        name.find("Contributions") == std::string::npos) {
      return name;
    }
  }
  return "";
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main.
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Expand either one ROOT file or an input directory into the files to process.
  const auto args = parse_args(argc, argv);
  std::vector<fs::path> input_files;
  const fs::path input_path(args.input_path);
  if (fs::is_regular_file(input_path) && input_path.extension() == ".root") {
    input_files.push_back(input_path);
  } else {
    input_files = rates::find_root_files(args.input_path);
  }
  if (input_files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_path << "\n";
    return 1;
  }

  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // This one decodes a cell ID into its corresponding position.
  const rates::HcalEndcapPInsertCellIDDecoder decoder;

  // This one maps a cell to a LFHCal-style channel/chip.
  const rates::InsertToLFHCALMapper mapper;

  // Keep one accumulator for each longitudinal layout and inner-ring treatment.
  rates::insert_analysis::physical_layers::Outputs physical_layers_with_inner_ring;
  rates::insert_analysis::physical_layers::Outputs physical_layers_no_inner_ring;
  rates::insert_analysis::lfhcal_segments::Outputs lfhcal_segments_with_inner_ring;
  rates::insert_analysis::lfhcal_segments::Outputs lfhcal_segments_no_inner_ring;

  rates::insert_analysis::physical_layers::with_inner_ring::init_outputs(physical_layers_with_inner_ring);
  rates::insert_analysis::physical_layers::no_inner_ring::init_outputs(physical_layers_no_inner_ring);
  rates::insert_analysis::lfhcal_segments::with_inner_ring::init_outputs(lfhcal_segments_with_inner_ring);
  rates::insert_analysis::lfhcal_segments::no_inner_ring::init_outputs(lfhcal_segments_no_inner_ring);

  std::string insert_collection_name;
  std::uint64_t n_events = 0;
  rates::FileProgress progress(input_files.size(), std::cerr);

  // Loop over all input files.
  for (const auto& path : input_files) {
    progress.tick();

    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Loop over all events in this file.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto event_data = reader.readEntry("events", event_index);
      if (!event_data) continue;

      podio::Frame frame(std::move(event_data));
      // Find the insert collection once, then require it in every later event.
      if (insert_collection_name.empty()) {
        insert_collection_name = find_insert_collection(frame);
      }
      if (insert_collection_name.empty() || !rates::has_collection(frame, insert_collection_name)) continue;
      ++n_events; // Count events that contribute to the rate normalization.

      // Copy each accepted hit into the compact form shared by the output layouts.
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(insert_collection_name);
      std::vector<rates::insert_analysis::EventHit> event_hits;
      event_hits.reserve(hits.size());

      // Loop over insert hits and discard cells outside the expected geometry.
      for (const auto& hit : hits) {
        const auto cell = decoder.cell(static_cast<std::uint64_t>(hit.getCellID()));
        if (cell.layer < 1 || cell.layer > rates::insert_analysis::physical_layers::kNLayers) continue;
        if (cell.side < 0 || cell.side > 1) continue;

        const auto position = hit.getPosition();
        event_hits.push_back(rates::insert_analysis::EventHit{
            cell.layer,
            cell.side,
            static_cast<double>(position.x),
            static_cast<double>(position.y),
            hit.getEnergy(),
        });
      }

      // Fold this event into each requested layout after every hit has been decoded.
      rates::insert_analysis::physical_layers::with_inner_ring::accumulate_event(
          physical_layers_with_inner_ring,
          event_hits,
          mapper);
      rates::insert_analysis::physical_layers::no_inner_ring::accumulate_event(
          physical_layers_no_inner_ring,
          event_hits,
          mapper);
      rates::insert_analysis::lfhcal_segments::with_inner_ring::accumulate_event(
          lfhcal_segments_with_inner_ring,
          event_hits,
          mapper);
      rates::insert_analysis::lfhcal_segments::no_inner_ring::accumulate_event(
          lfhcal_segments_no_inner_ring,
          event_hits,
          mapper);
    }
  }
  std::cerr << "\n";

  // Write every layout into the same ROOT file after the full sample is accumulated.
  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  rates::insert_analysis::physical_layers::with_inner_ring::write_output(
      output,
      physical_layers_with_inner_ring,
      n_events);
  rates::insert_analysis::physical_layers::no_inner_ring::write_output(
      output,
      physical_layers_no_inner_ring,
      n_events);
  rates::insert_analysis::lfhcal_segments::with_inner_ring::write_output(
      output,
      lfhcal_segments_with_inner_ring,
      n_events);
  rates::insert_analysis::lfhcal_segments::no_inner_ring::write_output(
      output,
      lfhcal_segments_no_inner_ring,
      n_events);

  output.Close();
  return 0;
}
