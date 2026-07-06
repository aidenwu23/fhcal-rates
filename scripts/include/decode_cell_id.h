#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string_view>

namespace rates {

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

struct LFHCALCellPosition {
  double x_mm = 0.0;
  double y_mm = 0.0;
};

struct LFHCALChipID {
  int moduleIDx = 0;
  int moduleIDy = 0;
  int module_type = 0;
  int module_half = 0;

  bool operator==(const LFHCALChipID& other) const {
    return moduleIDx == other.moduleIDx &&
           moduleIDy == other.moduleIDy &&
           module_type == other.module_type &&
           module_half == other.module_half;
  }
};

struct LFHCALChipIDHash {
  std::size_t operator()(const LFHCALChipID& chip) const;
};

class LFHCALCellIDDecoder {
 public:
  int get(std::uint64_t cell_id, std::string_view field) const;
  LFHCALChannelID channel(std::uint64_t cell_id) const;
  // Decode the chip (one virtual 4M module in the geometry view) seen by the hardware readout chip.
  LFHCALChipID decode_chip(std::uint64_t cell_id) const;
  bool is_passive(std::uint64_t cell_id) const;
  LFHCALCellPosition position(std::uint64_t cell_id) const;
  LFHCALCellPosition position(const LFHCALChipID& chip) const;
};

using LFHCALDecoder = LFHCALCellIDDecoder;

}  // namespace rates
