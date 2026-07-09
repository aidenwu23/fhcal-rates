#include "insert_to_lfhcal.h"

#include <cmath>

namespace rates {
namespace {

constexpr double kVirtualCellSizeMm = 50.0;

// Bin an x or y coordinate into the corresponding 5 cm wide virtual-channel bin.
int virtual_index(double coordinate_mm) {
  return static_cast<int>(std::floor(coordinate_mm / kVirtualCellSizeMm));
}

}  // namespace

std::size_t VirtualLFHCALChannelIDHash::operator()(const VirtualLFHCALChannelID& channel) const {
  std::size_t h = 0;
  const auto mix = [&](int value) {
    h ^= std::hash<int>{}(value) + 0x9e3779b9 + (h << 6) + (h >> 2);
  };

  mix(channel.layer);
  mix(channel.ix);
  mix(channel.iy);
  return h;
}

std::size_t VirtualLFHCALChipIDHash::operator()(const VirtualLFHCALChipID& chip) const {
  std::size_t h = 0;
  const auto mix = [&](int value) {
    h ^= std::hash<int>{}(value) + 0x9e3779b9 + (h << 6) + (h >> 2);
  };

  mix(chip.layer);
  mix(chip.ix);
  mix(chip.iy);
  return h;
}

VirtualLFHCALChannelID InsertToLFHCALMapper::channel(int layer, double x_mm, double y_mm) const {
  return VirtualLFHCALChannelID{
      layer,
      virtual_index(x_mm),
      virtual_index(y_mm),
  };
}

VirtualLFHCALChannelID InsertToLFHCALMapper::channel(const HcalEndcapPInsertCellID& cell,
                                                     double x_mm,
                                                     double y_mm) const {
  return channel(cell.layer, x_mm, y_mm);
}

VirtualLFHCALChipID InsertToLFHCALMapper::chip(const VirtualLFHCALChannelID& channel) const {
  return VirtualLFHCALChipID{
      1,
      // Integer division is used here, so 0 / 2 = 0 and 1 / 2 = 0 after truncation.
      // That makes neighboring channel columns share one chip column:
      // channel ix = 0,1 -> chip ix = 0; 2,3 -> 1; 4,5 -> 2; ...
      channel.ix / 2,
      // Do the same integer-division grouping in y so one chip covers a 2 by 2 block of channels.
      channel.iy / 2,
  };
}

}  // namespace rates
