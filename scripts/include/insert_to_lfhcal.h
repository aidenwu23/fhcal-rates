#pragma once

#include "decode_cell_id.h"

#include <cstddef>
#include <cstdint>
#include <functional>

namespace rates {

struct VirtualLFHCALChannelID {
  int layer = 0;
  int ix = 0;
  int iy = 0;

  bool operator==(const VirtualLFHCALChannelID& other) const {
    return layer == other.layer &&
           ix == other.ix &&
           iy == other.iy;
  }
};

struct VirtualLFHCALChannelIDHash {
  std::size_t operator()(const VirtualLFHCALChannelID& channel) const;
};

struct VirtualLFHCALChipID {
  int layer = 0;
  int ix = 0;
  int iy = 0;

  bool operator==(const VirtualLFHCALChipID& other) const {
    return layer == other.layer &&
           ix == other.ix &&
           iy == other.iy;
  }
};

struct VirtualLFHCALChipIDHash {
  std::size_t operator()(const VirtualLFHCALChipID& chip) const;
};

class InsertToLFHCALMapper {
 public:
  VirtualLFHCALChannelID channel(int layer, double x_mm, double y_mm) const;
  VirtualLFHCALChannelID channel(const HcalEndcapPInsertCellID& cell, double x_mm, double y_mm) const;
  VirtualLFHCALChipID chip(const VirtualLFHCALChannelID& channel) const;
};

}  // namespace rates
