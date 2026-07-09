/*

./build/insert_occupancy -i data/bkg_apr -o insert_plots/occupancy.root

*/

#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "channel.h"
#include "chip.h"
#include "decode_cell_id.h"
#include "mip.h"
#include "radius.h"
#include "side.h"
#include "insert_to_lfhcal.h"
#include "shared.h"
#include "utils.h"
#include "xy.h"

#include <array>
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

constexpr const char* kCollectionMatch = "HcalEndcapPInsert";

struct Args {
  std::string input_path;
  std::string output_file;
};

// ----------------------------------------------------------------------------------
// CLI
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

}  // namespace

// ----------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse args and expand the input into a file list.
  const auto args = parse_args(argc, argv);
  std::vector<fs::path> files;
  const fs::path input_path(args.input_path);
  if (fs::is_regular_file(input_path) && input_path.extension() == ".root") {
    files.push_back(input_path);
  } else {
    files = rates::find_root_files(args.input_path);
  }
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_path << "\n";
    return 1;
  }

  // Ensure the output directory exists before any event work starts.
  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Decodes the relevant bitfields out of the cell ID.
  const rates::HcalEndcapPInsertCellIDDecoder decoder;

  // Maps a cell onto either a virtual channel or virtual chip (like the lfhcal).
  // Takes the layer, x_mm, and y_mm.
  const rates::InsertToLFHCALMapper mapper;

  std::vector<rates::insert_occupancy::LayerAccum> layers;
  std::vector<rates::insert_occupancy::chip::LayerAccum> chip_layers;
  std::array<rates::insert_occupancy::mip::ChannelThresholdAccum, 16> channel_mip_products;
  std::array<rates::insert_occupancy::mip::ChipThresholdAccum, 16> chip_mip_products;
  std::array<rates::insert_occupancy::mip::DataThresholdAccum, 16> data_mip_products;
  std::array<rates::insert_occupancy::radius::ThresholdAccum, 3> radius_products;
  std::array<std::array<rates::insert_occupancy::side::ThresholdAccum, 16>, 2> side_products;
  rates::insert_occupancy::init_layer_accumulations(layers);
  rates::insert_occupancy::chip::init_layer_accumulations(chip_layers);

  std::string collection_name;
  std::uint64_t n_events = 0;
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
      ++n_events;

      // These maps hold the event-level virtual readout before it is folded into totals.
      std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash> event_hits;
      rates::insert_occupancy::mip::EventEnergyMap event_energy;
      std::array<rates::insert_occupancy::mip::EventEnergyMap, 2> side_event_energy;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(collection_name);

      // Loop over all hits.
      for (const auto& hit : hits) {

        // Decode cell.
        const auto cell = decoder.cell(static_cast<std::uint64_t>(hit.getCellID()));
        if (cell.layer < 1 || cell.layer > rates::insert_occupancy::kNLayers) continue;

        // Map the physical layer into one of the seven LFHCAL-like longitudinal segments.
        if (cell.side < 0 || cell.side > 1) continue;

        const int segment = rates::insert_occupancy::segment_index(cell.layer);
        if (segment < 0) continue;

        const auto position = hit.getPosition();

        // Map the hit into the virtual channel grid used for all later rate products.
        auto channel = mapper.channel(cell.layer, position.x, position.y);
        channel.layer = rates::insert_occupancy::kSegmentFirstLayers[segment];
        ++event_hits[channel];              // Count how many hits land in this virtual channel this event.
        event_energy[channel] += hit.getEnergy();  // Sum energy so threshold scans can be applied after the hit loop.
        side_event_energy[cell.side][channel] += hit.getEnergy();
      }

      // Convert the event-level virtual channels into each output product family.
      rates::insert_occupancy::xy::accumulate_event(layers, event_hits);
      rates::insert_occupancy::chip::accumulate_event(chip_layers, event_hits, mapper);
      rates::insert_occupancy::mip::accumulate_event(channel_mip_products, chip_mip_products, data_mip_products, event_energy, mapper);
      rates::insert_occupancy::radius::accumulate_event(radius_products, event_energy, mapper);
      rates::insert_occupancy::side::accumulate_event(side_products, side_event_energy, mapper);
      rates::insert_occupancy::fill_event_histograms(layers, event_hits);
    }
  }
  std::cerr << "\n";

  if (collection_name.empty()) {
    std::cerr << "Failed to find an insert hit collection\n";
    return 1;
  }
  if (n_events == 0) {
    std::cerr << "No events with insert hits found\n";
    return 1;
  }

  // Open the output file and write each product family into its own directory tree.
  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  rates::insert_occupancy::xy::write_output(output, layers, n_events);
  rates::insert_occupancy::channel::write_output(output, layers, n_events);
  rates::insert_occupancy::chip::write_output(output, chip_layers, n_events);
  rates::insert_occupancy::mip::write_output(output, channel_mip_products, chip_mip_products, data_mip_products, n_events);
  rates::insert_occupancy::radius::write_output(output, radius_products, n_events);
  rates::insert_occupancy::side::write_output(output, side_products, n_events);

  output.Close();
  return 0;
}
