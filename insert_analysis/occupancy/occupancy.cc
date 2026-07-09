/*

./build/insert_occupancy -i data/bkg_apr -o insert_plots/occupancy.root

*/

#include <TFile.h>
#include <TH1.h>

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

constexpr const char* kInsertCollectionMatch = "HcalEndcapPInsert";

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
    if (name.find(kInsertCollectionMatch) != std::string::npos &&
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

  // Ensure the output directory exists before any event work starts.
  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Decodes the relevant bitfields out of the cell ID.
  const rates::HcalEndcapPInsertCellIDDecoder decoder;

  // Maps a cell onto either a virtual channel or virtual chip (like the lfhcal).
  // Takes the layer, x_mm, and y_mm.
  const rates::InsertToLFHCALMapper mapper;

  std::vector<rates::insert_occupancy::SegmentSum> segment_sums;
  rates::insert_occupancy::chip::ChipSum chip_sum;
  std::array<rates::insert_occupancy::mip::ChannelThresholdSum, 16> channel_threshold_sums;
  std::array<rates::insert_occupancy::mip::ChipThresholdSum, 16> chip_threshold_sums;
  std::array<rates::insert_occupancy::mip::DataThresholdSum, 16> data_threshold_sums;
  std::array<rates::insert_occupancy::radius::ThresholdSum, 3> radius_threshold_sums;
  std::array<std::array<rates::insert_occupancy::side::ThresholdSum, 16>, 2> side_threshold_sums;
  rates::insert_occupancy::init_segment_sums(segment_sums);
  rates::insert_occupancy::chip::init_chip_sum(chip_sum);

  std::string insert_collection_name;
  std::uint64_t n_events = 0;
  rates::FileProgress progress(input_files.size(), std::cerr);

  // Loop over all files in the input.
  for (const auto& path : input_files) {
    progress.tick();

    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Per file, loop over all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto event_data = reader.readEntry("events", event_index);
      if (!event_data) continue;

      podio::Frame frame(std::move(event_data));

      if (insert_collection_name.empty()) {
        insert_collection_name = find_insert_collection(frame);
      }
      if (insert_collection_name.empty() || !rates::has_collection(frame, insert_collection_name)) continue;
      ++n_events;

      // These maps hold the event-level virtual readout before it is folded into totals.
      std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash> channel_hit_counts;
      rates::insert_occupancy::mip::EventEnergyMap channel_energy_sum;
      std::array<rates::insert_occupancy::mip::EventEnergyMap, 2> side_channel_energy_sum;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(insert_collection_name);

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
        ++channel_hit_counts[channel];
        channel_energy_sum[channel] += hit.getEnergy();
        side_channel_energy_sum[cell.side][channel] += hit.getEnergy();
      }

      // Build the fixed 0.5-MIP virtual-channel view used by the non-scan occupancy products.
      std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash> channel_hit_counts_above_0p5_mip;
      for (const auto& [channel, energy] : channel_energy_sum) {
        const int segment = rates::insert_occupancy::segment_index(channel.layer);
        if (energy <= 0.5 * rates::insert_occupancy::channel_mip_energy_gev(segment)) continue;
        auto it = channel_hit_counts.find(channel);
        if (it == channel_hit_counts.end()) continue;
        channel_hit_counts_above_0p5_mip[channel] = it->second;
      }

      // Convert the event-level virtual channels into each output product family.
      rates::insert_occupancy::xy::accumulate_event(segment_sums, channel_hit_counts_above_0p5_mip);
      rates::insert_occupancy::chip::accumulate_event(chip_sum, channel_hit_counts_above_0p5_mip, mapper);
      rates::insert_occupancy::mip::accumulate_event(channel_threshold_sums, chip_threshold_sums, data_threshold_sums, channel_energy_sum, mapper);
      rates::insert_occupancy::radius::accumulate_event(radius_threshold_sums, channel_energy_sum, mapper);
      rates::insert_occupancy::side::accumulate_event(side_threshold_sums, side_channel_energy_sum, mapper);
      rates::insert_occupancy::fill_event_histograms(segment_sums, channel_hit_counts_above_0p5_mip);
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

  // Open the output file and write each product family into its own directory tree.
  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  rates::insert_occupancy::xy::write_output(output, segment_sums, n_events);
  rates::insert_occupancy::channel::write_output(output, segment_sums, n_events);
  rates::insert_occupancy::chip::write_output(output, chip_sum, n_events);
  rates::insert_occupancy::mip::write_output(output, channel_threshold_sums, chip_threshold_sums, data_threshold_sums, n_events);
  rates::insert_occupancy::radius::write_output(output, radius_threshold_sums, n_events);
  rates::insert_occupancy::side::write_output(output, side_threshold_sums, n_events);

  output.Close();
  return 0;
}
