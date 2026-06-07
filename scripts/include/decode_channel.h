#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string_view>

namespace br {

struct LFHCALChannelID {
  int moduleIDx = 0;
  int moduleIDy = 0;
  int moduletype = 0;
  int towerx = 0;
  int towery = 0;
  int rlayerz = 0;

  bool operator==(const LFHCALChannelID& other) const {
    return moduleIDx == other.moduleIDx &&
           moduleIDy == other.moduleIDy &&
           moduletype == other.moduletype &&
           towerx == other.towerx &&
           towery == other.towery &&
           rlayerz == other.rlayerz;
  }
};

struct LFHCALChannelIDHash {
  std::size_t operator()(const LFHCALChannelID& channel) const;
};

class LFHCALDecoder {
 public:
  int get(std::uint64_t cell_id, std::string_view field) const;
  LFHCALChannelID channel(std::uint64_t cell_id) const;
};

}  // namespace br
