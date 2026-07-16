#include "insert_to_lfhcal.h"

#include <cmath>

namespace rates {
namespace {

constexpr double kVirtualCellSizeMm = 50.0;
constexpr double kPizzaCenterXMM = 0.5 * (-3.95 - 21.5);
constexpr double kPi = 3.14159265358979323846;
constexpr double kLeftBoundary0 = 0.9127426518583035;
constexpr double kLeftBoundary1 = 0.5 * kPi;
constexpr double kLeftBoundary2 = 2.2288500017314896;
constexpr double kRightBoundary0 = 0.8133195162331286;
constexpr double kRightBoundary1 = 2.3282731373566645;

// Bin an x or y coordinate into the corresponding 5 cm wide virtual-channel bin.
int virtual_index(double coordinate_mm) {
  return static_cast<int>(std::floor(coordinate_mm / kVirtualCellSizeMm));
}

// Group pairs of channel indices while preserving the pair boundary across zero.
int chip_index(int channel_index) {
  return channel_index >= 0 ? channel_index / 2 : (channel_index - 1) / 2;
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

std::size_t VirtualLFHCALPizzaChipIDHash::operator()(const VirtualLFHCALPizzaChipID& chip) const {
  std::size_t h = 0;
  const auto mix = [&](int value) {
    h ^= std::hash<int>{}(value) + 0x9e3779b9 + (h << 6) + (h >> 2);
  };

  mix(chip.layer);
  mix(chip.side);
  mix(chip.region);
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
      // Neighboring channel columns share one chip column on both sides of zero.
      chip_index(channel.ix),
      // Do the same integer-division grouping in y so one chip covers a 2 by 2 block of channels.
      chip_index(channel.iy),
  };
}

VirtualLFHCALPizzaChipID InsertToLFHCALMapper::pizza_chip(const VirtualLFHCALChannelID& channel,
                                                          int side) const {

  // Convert channel position into an angular value.
  const double x_mm = (static_cast<double>(channel.ix) + 0.5) * kVirtualCellSizeMm;
  const double y_mm = (static_cast<double>(channel.iy) + 0.5) * kVirtualCellSizeMm;
  const double angle = std::atan2(y_mm, x_mm - kPizzaCenterXMM);

  // Give each side its own angle from 0 to pi:
  // side 0 runs from the top edge, around the left, to the bottom edge;
  // side 1 runs from the bottom edge, around the right, to the top edge.
  double side_angle = side == 0 ? angle - 0.5 * kPi : angle + 0.5 * kPi;

  // Move negative results into the equivalent 0 to 2*pi range.
  if (side_angle < 0.0) side_angle += 2.0 * kPi;

  // A channel slightly beyond its side's half-circle is assigned to the nearest end of that side.
  if (side_angle > 1.5 * kPi) side_angle = 0.0;
  else if (side_angle > kPi) side_angle = kPi;

  int region = 0;
  if (side == 0) {
    if (side_angle >= 0.0 && side_angle < kLeftBoundary0) region = 0;
    else if (side_angle >= kLeftBoundary0 && side_angle < kLeftBoundary1) region = 1;
    else if (side_angle >= kLeftBoundary1 && side_angle < kLeftBoundary2) region = 2;
    else if (side_angle >= kLeftBoundary2 && side_angle <= kPi) region = 3;
  } else {
    if (side_angle >= 0.0 && side_angle < kRightBoundary0) region = 0;
    else if (side_angle >= kRightBoundary0 && side_angle < kRightBoundary1) region = 1;
    else if (side_angle >= kRightBoundary1 && side_angle <= kPi) region = 2;
  }

  return VirtualLFHCALPizzaChipID{
      channel.layer,
      side,
      region,
  };
}

}  // namespace rates
