/*

./build/lfhcal_occupancy -i data/bkg_july -o lfhcal_plots/occupancy.root

*/

#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>
#include <edm4hep/MCParticle.h>

#include "decode_cell_id.h"
#include "classify_hit.h"
#include "event_channel.h"
#include "occupancy/include/family.h"
#include "utils.h"

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kHitCollection = "LFHCALHits";

struct Args {
  std::string input_path;
  std::string output_file;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT.root|DIR -o OUTPUT.root\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;
  // Read the input and output paths from the command line.
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

}  // namespace

int main(int argc, char* argv[]) {
  // Configure ROOT and collect the input files.
  TH1::AddDirectory(false);
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetTitleW(0.8);

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

  // Prepare the output path and full-sample accumulators.
  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  const rates::LFHCALCellIDDecoder decoder;
  rates::lfhcal_analysis::occupancy::Outputs outputs;
  rates::lfhcal_analysis::occupancy::init_outputs(outputs);
  std::uint64_t n_events = 0;
  rates::FileProgress progress(input_files.size(), std::cerr);

  // Loop over input files.
  for (const auto& path : input_files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Loop over events in this file.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto event_data = reader.readEntry("events", event_index);
      if (!event_data) continue;
      podio::Frame frame(std::move(event_data));
      if (!rates::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      // Group this event's active hits by readout channel.
      std::unordered_map<rates::LFHCALChannelID,
                         rates::lfhcal_analysis::EventChannel,
                         rates::LFHCALChannelIDHash> event_channel_map;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop over all hits.
      for (const auto& hit : hits) {

        // Decode hit cell ID.
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;
        const auto channel_id = decoder.channel(cell_id);
        if (channel_id.rlayerz < 0 || channel_id.rlayerz >= rates::lfhcal_analysis::occupancy::kNReadoutLayers) continue;

        auto& event_channel = event_channel_map[channel_id];
        // Initialize channel metadata when the channel first appears.
        if (event_channel.hit_count == 0) {
          const auto position = decoder.position(cell_id);
          event_channel.channel = channel_id;
          event_channel.chip = decoder.decode_chip(cell_id);
          event_channel.x_mm = position.x_mm;
          event_channel.y_mm = position.y_mm;
        }
        event_channel.energy_gev += hit.getEnergy();
        
        // Sum the tracked background contributions within this channel.
        for (const auto& contribution : hit.getContributions()) {
          if (contribution.getEnergy() <= 0.0) continue;
          const int origin = rates::origin_index(contribution.getParticle().getGeneratorStatus());
          if (origin >= 0 && origin < static_cast<int>(event_channel.energy_by_origin.size())) {
            event_channel.energy_by_origin[origin] += contribution.getEnergy();
          }
        }
        ++event_channel.hit_count;
      }

      // Convert the channel map and send the completed event to every analysis family.
      rates::lfhcal_analysis::occupancy::EventChannels event_channels;
      event_channels.reserve(event_channel_map.size());
      for (const auto& [channel_id, event_channel] : event_channel_map) {
        (void)channel_id;
        event_channels.push_back(event_channel);
      }
      rates::lfhcal_analysis::occupancy::accumulate_event(outputs, event_channels, decoder);
    }
  }
  std::cerr << "\n";

  // Validate the event count and write the consolidated output file.
  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }
  rates::lfhcal_analysis::occupancy::write_output(output, outputs, decoder, n_events);
  output.Close();
  return 0;
}
