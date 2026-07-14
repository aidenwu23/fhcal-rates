/*

root -l -b -q 'insert_analysis/count_channels.C("data/bkg_apr")'

*/

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

void count_channels(const char* input_dir = "data/bkg_apr") {
  std::vector<std::string> files;
  for (const auto& entry : fs::recursive_directory_iterator(input_dir)) {
    if (entry.is_regular_file() && entry.path().extension() == ".root") {
      files.push_back(entry.path().string());
    }
  }
  std::sort(files.begin(), files.end());

  std::array<std::array<std::unordered_set<std::uint64_t>, 2>, 60> channels_by_layer_and_side;
  std::string collection_name;

  for (const auto& file : files) {
    podio::ROOTReader reader;
    reader.openFile(file);

    for (std::size_t event = 0; event < reader.getEntries("events"); ++event) {
      auto data = reader.readEntry("events", event);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (collection_name.empty()) {
        for (const auto& name : frame.getAvailableCollections()) {
          if (name.find("HcalEndcapPInsert") != std::string::npos &&
              name.find("Contributions") == std::string::npos) {
            collection_name = name;
            break;
          }
        }
      }
      if (collection_name.empty()) continue;

      const auto& collections = frame.getAvailableCollections();
      if (std::find(collections.begin(), collections.end(), collection_name) == collections.end()) continue;

      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(collection_name);
      for (const auto& hit : hits) {
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        const int side = static_cast<int>((cell_id >> 8) & 0x1);
        const int layer = static_cast<int>((cell_id >> 9) & 0xff);
        if (layer < 1 || layer > 60) continue;

        const auto position = hit.getPosition();
        const int ix = static_cast<int>(std::floor(static_cast<double>(position.x) / 50.0));
        const int iy = static_cast<int>(std::floor(static_cast<double>(position.y) / 50.0));
        const std::uint64_t channel =
            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(ix)) << 32) |
            static_cast<std::uint32_t>(iy);
        channels_by_layer_and_side[layer - 1][side].insert(channel);
      }
    }
  }

  std::cout << std::setw(5) << "layer"
            << std::setw(10) << "side 0"
            << std::setw(10) << "side 1" << "\n";
  for (int layer = 1; layer <= 60; ++layer) {
    std::cout << std::setw(5) << layer
              << std::setw(10) << channels_by_layer_and_side[layer - 1][0].size()
              << std::setw(10) << channels_by_layer_and_side[layer - 1][1].size() << "\n";
  }
}
