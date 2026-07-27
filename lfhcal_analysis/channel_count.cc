/*

./build/channel_count -i data/bkg_july

*/

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "utils.h"

#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_set>

namespace {

constexpr const char* kHitCollection = "LFHCALHits";
constexpr int kNLayers = 7;

struct Args {
  std::string input_dir;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT_DIR\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_dir = argv[++i];
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input_dir.empty()) {
    usage(argv[0]);
    std::exit(1);
  }

  return args;
}

}  // namespace

int main(int argc, char* argv[]) {
  const auto args = parse_args(argc, argv);
  const auto files = rates::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  const rates::LFHCALCellIDDecoder decoder;
  std::unordered_set<rates::LFHCALChannelID, rates::LFHCALChannelIDHash> all_channels;
  rates::FileProgress progress(files.size(), std::cerr);

  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!rates::has_collection(frame, kHitCollection)) continue;

      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        const auto channel = decoder.channel(cell_id);
        if (channel.rlayerz < 0 || channel.rlayerz >= kNLayers) continue;

        all_channels.insert(channel);
      }
    }
  }

  std::cout << "\n";
  std::cout << all_channels.size() << "\n";

  return 0;
}
