#include "decode_cell_id.h"

#include <stdexcept>

namespace br {
namespace {

/*
system:8,moduleIDx:6,moduleIDy:6,moduletype:1,passive:1,
towerx:2,towery:1,rlayerz:4,layerz:4
*/

// system
constexpr int kSystemOffset = 0;
constexpr int kSystemWidth = 8;

// moduleIDx
constexpr int kModuleIDxOffset = kSystemOffset + kSystemWidth;
constexpr int kModuleIDxWidth = 6;

// moduleIDy
constexpr int kModuleIDyOffset = kModuleIDxOffset + kModuleIDxWidth;
constexpr int kModuleIDyWidth = 6;

// moduletype (0 = 8M module, 1 = 4M module).
constexpr int kModuletypeOffset = kModuleIDyOffset + kModuleIDyWidth;
constexpr int kModuletypeWidth = 1;

// passive
constexpr int kPassiveOffset = kModuletypeOffset + kModuletypeWidth;
constexpr int kPassiveWidth = 1;

// towerx
constexpr int kTowerxOffset = kPassiveOffset + kPassiveWidth;
constexpr int kTowerxWidth = 2;

// towery
constexpr int kToweryOffset = kTowerxOffset + kTowerxWidth;
constexpr int kToweryWidth = 1;

// rlayerz
constexpr int kRlayerzOffset = kToweryOffset + kToweryWidth;
constexpr int kRlayerzWidth = 4;

// layerz
constexpr int kLayerzOffset = kRlayerzOffset + kRlayerzWidth;
constexpr int kLayerzWidth = 4;

// Builld a mask with "width" low bits set.
std::uint64_t bit_mask(int width) {
  return (1ULL << width) - 1ULL;
}

// Decode one DD4hep bitfield:
//   1. shift the target field down to bit 0
//   2. mask away all other bits
//   3. cast to int for downstream indexing
int decode_bits(std::uint64_t cell_id, int offset, int width) {
  return static_cast<int>((cell_id >> offset) & bit_mask(width));
}

/*
Derived plotting lattice for 8M x.

LFHCAL_geo.cpp uses:
  moduleIDx = (pos8M.x + 270) / 10
  placement x = -pos8M.x

Therefore:
  module center x_cm = 270 - 10 * moduleIDx

Ideal 8M tower centers for x are:
  +7.5, +2.5, -2.5, -7.5 cm
which is:
  7.5 - 5 * tower_x

Divide by 2.5 cm to get an integer representative of the x coordinate (lattice):
  (270 - 10*m + 7.5 - 5*t) / 2.5
  = 108 - 4*m + 3 - 2*t
  = 111 - 4*m - 2*t
*/
int x_key_8m(int module_id_x, int tower_x) {
  return 111 - 4 * module_id_x - 2 * tower_x;
}

/*
Derived plotting lattice for 4M x.

LFHCAL_geo.cpp uses:
  moduleIDx = (pos4M.x + 265) / 10
  placement x = -pos4M.x

Therefore:
  module center x_cm = 265 - 10 * moduleIDx

Ideal 4M tower centers are:
  +2.5, -2.5 cm
which is:
  2.5 - 5 * tower_x

Divide by 2.5 cm:
  (265 - 10*m + 2.5 - 5*t) / 2.5
  = 106 - 4*m + 1 - 2*t
  = 107 - 4*m - 2*t
*/
int x_key_4m(int module_id_x, int tower_x) {
  return 107 - 4 * module_id_x - 2 * tower_x;
}

/*
Derived plotting lattice for y.

LFHCAL_geo.cpp uses:
  moduleIDy = (pos.y + 265) / 10
  placement y = -pos.y

Therefore:
  module center y_cm = 265 - 10 * moduleIDy

Ideal y tower centers are:
  +2.5, -2.5 cm
which is:
  2.5 - 5 * tower_y

Divide by 2.5 cm:
  (265 - 10*m + 2.5 - 5*t) / 2.5
  = 107 - 4*m - 2*t
*/
int y_key(int module_id_y, int tower_y) {
  return 107 - 4 * module_id_y - 2 * tower_y;
}

}  // namespace


// Combine the channel fields into one hash value.
// This lets LFHCALChannelID work as a lookup value in unordered_map / unordered_set.
std::size_t LFHCALChannelIDHash::operator()(const LFHCALChannelID& channel) const {

  std::size_t h = 0;

  // Hash-combine pattern:
  //   hash current field
  //   add a mixing constant
  //   mix with previous accumulated hash
  //
  // The constant 0x9e3779b9 is commonly used in boost-style hash combining.  
  const auto mix = [&](int value) {
    h ^= std::hash<int>{}(value) + 0x9e3779b9 + (h << 6) + (h >> 2);
  };

  // Channel identity = transverse cell location & readout layer. (via rlayerz).
  mix(channel.moduleIDx);
  mix(channel.moduleIDy);
  mix(channel.moduletype);
  mix(channel.towerx);
  mix(channel.towery);
  mix(channel.rlayerz);
  return h;
}

// Decode a LFHCAL DD4hep field from the packed 64-bit cell_id.
// Field offsets and widths come from the LFHCAL readout string.
int LFHCALCellIDDecoder::get(std::uint64_t cell_id, std::string_view field) const {
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
LFHCALChannelID LFHCALCellIDDecoder::channel(std::uint64_t cell_id) const {
  return LFHCALChannelID{
      get(cell_id, "moduleIDx"),
      get(cell_id, "moduleIDy"),
      get(cell_id, "moduletype"),
      get(cell_id, "towerx"),
      get(cell_id, "towery"),
      get(cell_id, "rlayerz"),
  };
}

bool LFHCALCellIDDecoder::is_passive(std::uint64_t cell_id) const {
  return get(cell_id, "passive") != 0;
}

/*
Convert a cell_id into a transverse plotting position.

The x/y lattice values above are in units of 2.5 cm.
Multiplying by 25.0 converts those units to mm:

  1 lattice unit = 2.5 cm = 25 mm

If LFHCALCellPosition expects cm, use 2.5 here.
If LFHCALCellPosition expects mm, keep 25.0.
*/
LFHCALCellPosition LFHCALCellIDDecoder::position(std::uint64_t cell_id) const {
  const int module_id_x = get(cell_id, "moduleIDx");
  const int module_id_y = get(cell_id, "moduleIDy");
  const int module_type = get(cell_id, "moduletype");
  const int tower_x = get(cell_id, "towerx");
  const int tower_y = get(cell_id, "towery");

  // moduletype selects the x-origin and tower layout:
  //   0 -> 8M: x origin 270 cm, four x towers
  //   1 -> 4M: x origin 265 cm, two x towers
  //
  // y uses the same two-row tower structure for both module families.

  const int x_key = module_type == 0 ? x_key_8m(module_id_x, tower_x) : x_key_4m(module_id_x, tower_x);
  const int y_lattice_key = y_key(module_id_y, tower_y);

  return LFHCALCellPosition{
      25.0 * static_cast<double>(x_key),
      25.0 * static_cast<double>(y_lattice_key),
  };
}

}  // namespace br
