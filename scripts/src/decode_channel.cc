#include "decode_channel.h"

#include <stdexcept>

namespace br {
namespace {

// System:starts at bit 0, 8 bits wide.
constexpr int kSystemOffset = 0;
constexpr int kSystemWidth = 8;

// moduleIDx: starts after system (at bit 8), 6 bits wide.
constexpr int kModuleIDxOffset = kSystemOffset + kSystemWidth;
constexpr int kModuleIDxWidth = 6;

// moduleIDy: starts after moduleIDx (bit 14), 6 bits wide.
constexpr int kModuleIDyOffset = kModuleIDxOffset + kModuleIDxWidth;
constexpr int kModuleIDyWidth = 6;

// moduletype: starts after moduleIDy (bit 20), 1 bit wide.
constexpr int kModuletypeOffset = kModuleIDyOffset + kModuleIDyWidth;
constexpr int kModuletypeWidth = 1;

// passive: starts after moduletype (bit 21), 1 bit wide.
constexpr int kPassiveOffset = kModuletypeOffset + kModuletypeWidth;
constexpr int kPassiveWidth = 1;

// towerx: starts after passsive (bit 22), 2 bits wide; tower index in x inside a module.
constexpr int kTowerxOffset = kPassiveOffset + kPassiveWidth;
constexpr int kTowerxWidth = 2;

// towery: starts after towerx (bit 24), 1 bit wide; tower index in y inside a module.
constexpr int kToweryOffset = kTowerxOffset + kTowerxWidth;
constexpr int kToweryWidth = 1;

// rlayerz: start after towery (bit 25), 4 bits wide: readout-layer index (total of 7).
constexpr int kRlayerzOffset = kToweryOffset + kToweryWidth;
constexpr int kRlayerzWidth = 4;

// layerz: start after rlayerz (bit 29), 4 bits wide; full geometry layer index.
constexpr int kLayerzOffset = kRlayerzOffset + kRlayerzWidth;
constexpr int kLayerzWidth = 4;

// Builds a mask of width 1-bits to extract a desired bitfield.
std::uint64_t bit_mask(int width) {
  return (1ULL << width) - 1ULL;
}

// Shift cell_id right by a certain offset so the wanted field lands at the far right, then apply the mask
// to extract that field, and convert the result to int.
int decode_bits(std::uint64_t cell_id, int offset, int width) {
  return static_cast<int>((cell_id >> offset) & bit_mask(width));
}

}  // namespace

// Build one hash number from the relevant channel fields to create a lookup key.
std::size_t LFHCALChannelIDHash::operator()(const LFHCALChannelID& channel) const {

  // Init the hash accummulator.
  std::size_t h = 0;

  // Create a function that can use outer varibles by reference.
  const auto mix = [&](int value) {

    // xor-assign a certain hash field into h and (likely) make a unique key by scrambling the bits more.
    h ^= std::hash<int>{}(value) + 0x9e3779b9 + (h << 6) + (h >> 2);
  };

  // Accumulate all relevant fields for a channel and return the hash value.
  mix(channel.moduleIDx);
  mix(channel.moduleIDy);
  mix(channel.moduletype);
  mix(channel.towerx);
  mix(channel.towery);
  mix(channel.rlayerz);
  return h;
}

// Take a cell_id and decode the relevant field out of it.
int LFHCALDecoder::get(std::uint64_t cell_id, std::string_view field) const {
  if (field == "system") return decode_bits(cell_id, kSystemOffset, kSystemWidth);
  if (field == "moduleIDx") return decode_bits(cell_id, kModuleIDxOffset, kModuleIDxWidth);
  if (field == "moduleIDy") return decode_bits(cell_id, kModuleIDyOffset, kModuleIDyWidth);
  if (field == "moduletype") return decode_bits(cell_id, kModuletypeOffset, kModuletypeWidth);
  if (field == "passive") return decode_bits(cell_id, kPassiveOffset, kPassiveWidth);
  if (field == "towerx") return decode_bits(cell_id, kTowerxOffset, kTowerxWidth);
  if (field == "towery") return decode_bits(cell_id, kToweryOffset, kToweryWidth);
  if (field == "rlayerz") return decode_bits(cell_id, kRlayerzOffset, kRlayerzWidth);
  if (field == "layerz") return decode_bits(cell_id, kLayerzOffset, kLayerzWidth);
  throw std::invalid_argument("Unknown LFHCAL cellID field");
}

// Assign a cell_id to a channel_id (all cells that share transverse location in a layer).
LFHCALChannelID LFHCALDecoder::channel(std::uint64_t cell_id) const {
  return LFHCALChannelID{
      get(cell_id, "moduleIDx"),
      get(cell_id, "moduleIDy"),
      get(cell_id, "moduletype"),
      get(cell_id, "towerx"),
      get(cell_id, "towery"),
      get(cell_id, "rlayerz"),
  };
}

}  // namespace br
