/*

root -l -b -q 'insert_analysis/count_channels.C("data/bkg_july")'

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
  constexpr double kVirtualCellSizeMM = 50.0;
  constexpr double kPizzaCenterXMM = 0.5 * (-3.95 - 21.5);
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kLeftBoundary0 = 0.9127426518583035;
  constexpr double kLeftBoundary1 = 0.5 * kPi;
  constexpr double kLeftBoundary2 = 2.2288500017314896;
  constexpr double kRightBoundary0 = 0.8133195162331286;
  constexpr double kRightBoundary1 = 2.3282731373566645;

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
        const int ix = static_cast<int>(std::floor(static_cast<double>(position.x) / kVirtualCellSizeMM));
        const int iy = static_cast<int>(std::floor(static_cast<double>(position.y) / kVirtualCellSizeMM));
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

  const auto channel_angle = [&](std::uint64_t channel) {
    const int ix = static_cast<std::int32_t>(channel >> 32);
    const int iy = static_cast<std::int32_t>(channel & 0xffffffff);
    const double x_mm = (static_cast<double>(ix) + 0.5) * kVirtualCellSizeMM;
    const double y_mm = (static_cast<double>(iy) + 0.5) * kVirtualCellSizeMM;
    return std::atan2(y_mm, x_mm - kPizzaCenterXMM);
  };

  const auto pizza_chip = [&](std::uint64_t channel, int side) {
    const double angle = channel_angle(channel);
    double side_angle = side == 0 ? angle - 0.5 * kPi : angle + 0.5 * kPi;
    if (side_angle < 0.0) side_angle += 2.0 * kPi;
    if (side_angle > 1.5 * kPi) side_angle = 0.0;
    else if (side_angle > kPi) side_angle = kPi;

    int chip = 0;
    if (side == 0) {
      if (side_angle >= 0.0 && side_angle < kLeftBoundary0) chip = 0;
      else if (side_angle >= kLeftBoundary0 && side_angle < kLeftBoundary1) chip = 1;
      else if (side_angle >= kLeftBoundary1 && side_angle < kLeftBoundary2) chip = 2;
      else if (side_angle >= kLeftBoundary2 && side_angle <= kPi) chip = 3;
    } else {
      if (side_angle >= 0.0 && side_angle < kRightBoundary0) chip = 0;
      else if (side_angle >= kRightBoundary0 && side_angle < kRightBoundary1) chip = 1;
      else if (side_angle >= kRightBoundary1 && side_angle <= kPi) chip = 2;
    }
    return chip;
  };

  std::cout << "\nApproximate pizza-chip channels\n";
  std::cout << std::setw(5) << "layer"
            << std::setw(8) << "left 0"
            << std::setw(8) << "left 1"
            << std::setw(8) << "left 2"
            << std::setw(8) << "left 3"
            << std::setw(9) << "right 0"
            << std::setw(9) << "right 1"
            << std::setw(9) << "right 2" << "\n";
  for (int layer = 0; layer < 60; ++layer) {
    std::array<std::array<int, 4>, 2> chip_channels{};
    for (int side = 0; side < 2; ++side) {
      for (const auto channel : channels_by_layer_and_side[layer][side]) {
        ++chip_channels[side][pizza_chip(channel, side)];
      }
    }

    std::cout << std::setw(5) << layer + 1;
    for (int chip = 0; chip < 4; ++chip) std::cout << std::setw(8) << chip_channels[0][chip];
    for (int chip = 0; chip < 3; ++chip) std::cout << std::setw(9) << chip_channels[1][chip];
    std::cout << "\n";
  }
}
