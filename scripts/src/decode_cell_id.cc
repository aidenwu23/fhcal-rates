#include "decode_cell_id.h"

#include <cstdint>
#include <stdexcept>

namespace rates {
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

int decode_signed_bits(std::uint64_t cell_id, int offset, int width) {
  const std::uint64_t raw = (cell_id >> offset) & bit_mask(width);
  const std::uint64_t sign_bit = 1ULL << (width - 1);
  if ((raw & sign_bit) == 0) {
    return static_cast<int>(raw);
  }

  const std::uint64_t extended = raw | ~bit_mask(width);
  return static_cast<int>(static_cast<std::int64_t>(extended));
}

/*
LFHCAL_geo.cpp uses:
  moduleIDx = (pos8M.x + 270) / 10                      (line 1122)
  placement x = -pos8M.x 

Therefore the 8M module center in world coordinates is:
  x_cm = 270 - 10 * moduleIDx

The four ideal 8M tower centers are:                    (lines 306-309)
  +7.5, +2.5, -2.5, -7.5 cm
which is:
  7.5 - 5 * towerx
*/
double x_cm_8m(int module_id_x, int tower_x) {
  return (270.0 - 10.0 * static_cast<double>(module_id_x)) + (7.5 - 5.0 * static_cast<double>(tower_x));
}

/*
LFHCAL_geo.cpp uses:
  moduleIDx = (pos4M.x + 265) / 10                      (line 1151)
  placement x = -pos4M.x

Therefore the 4M module center in world coordinates is:
  x_cm = 265 - 10 * moduleIDx

The two ideal 4M tower centers are:                     (lines 434-435)
  +2.5, -2.5 cm
which is:
  2.5 - 5 * towerx
*/
double x_cm_4m(int module_id_x, int tower_x) {
  return (265.0 - 10.0 * static_cast<double>(module_id_x)) + (2.5 - 5.0 * static_cast<double>(tower_x));
}

/*
LFHCAL_geo.cpp uses:
  moduleIDy = (pos.y + 265) / 10                        (lines 1123, 1152)
  placement y = -pos.y

Therefore the module center in world coordinates is:
  y_cm = 265 - 10 * moduleIDy

Both module families use two y towers with centers:     (lines 310-312, 436-437)
  +2.5, -2.5 cm
which is:
  2.5 - 5 * towery
*/
double y_cm(int module_id_y, int tower_y) {
  return (265.0 - 10.0 * static_cast<double>(module_id_y)) + (2.5 - 5.0 * static_cast<double>(tower_y));
}

double chip_x_cm(int module_id_x, int module_type, int module_half) {
  if (module_type == 0) {
    return (270.0 - 10.0 * static_cast<double>(module_id_x)) + (module_half == 0 ? 5.0 : -5.0);
  }
  return 265.0 - 10.0 * static_cast<double>(module_id_x);
}

double chip_y_cm(int module_id_y) {
  return 265.0 - 10.0 * static_cast<double>(module_id_y);
}

// HcalEndcapPInsert readout:
// system:8,side:1,layer:8,slice:7,x:32:-16,y:-16
//
// The insert geometry aligns x to bit 32, leaving 8 unused bits after slice.
constexpr int kInsertSystemOffset = 0;
constexpr int kInsertSystemWidth = 8;

constexpr int kInsertSideOffset = kInsertSystemOffset + kInsertSystemWidth;
constexpr int kInsertSideWidth = 1;

constexpr int kInsertLayerOffset = kInsertSideOffset + kInsertSideWidth;
constexpr int kInsertLayerWidth = 8;

constexpr int kInsertSliceOffset = kInsertLayerOffset + kInsertLayerWidth;
constexpr int kInsertSliceWidth = 7;

constexpr int kInsertXOffset = 32;
constexpr int kInsertXWidth = 16;

constexpr int kInsertYOffset = kInsertXOffset + kInsertXWidth;
constexpr int kInsertYWidth = 16;

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

std::size_t LFHCALChipIDHash::operator()(const LFHCALChipID& chip) const {
  std::size_t h = 0;
  const auto mix = [&](int value) {
    h ^= std::hash<int>{}(value) + 0x9e3779b9 + (h << 6) + (h >> 2);
  };

  mix(chip.moduleIDx);
  mix(chip.moduleIDy);
  mix(chip.module_type);
  mix(chip.module_half);
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

LFHCALChipID LFHCALCellIDDecoder::decode_chip(std::uint64_t cell_id) const {
  const int module_id_x = get(cell_id, "moduleIDx");
  const int module_id_y = get(cell_id, "moduleIDy");
  const int module_type = get(cell_id, "moduletype");
  const int tower_x = get(cell_id, "towerx");

  return LFHCALChipID{
      module_id_x,
      module_id_y,
      module_type,
      module_type == 0 ? (tower_x < 2 ? 0 : 1) : 0,
  };
}

bool LFHCALCellIDDecoder::is_passive(std::uint64_t cell_id) const {
  return get(cell_id, "passive") != 0;
}

LFHCALCellPosition LFHCALCellIDDecoder::position(std::uint64_t cell_id) const {
  const int module_id_x = get(cell_id, "moduleIDx");
  const int module_id_y = get(cell_id, "moduleIDy");
  const int module_type = get(cell_id, "moduletype");
  const int tower_x = get(cell_id, "towerx");
  const int tower_y = get(cell_id, "towery");

  // moduletype selects the x-origin and x-tower layout:
  //   0 -> 8M: x origin 270 cm, four x towers
  //   1 -> 4M: x origin 265 cm, two x towers
  const double x_cm = module_type == 0 ? x_cm_8m(module_id_x, tower_x) : x_cm_4m(module_id_x, tower_x);
  const double y_cm_value = y_cm(module_id_y, tower_y);

  return LFHCALCellPosition{
      10.0 * x_cm,
      10.0 * y_cm_value,
  };
}

LFHCALCellPosition LFHCALCellIDDecoder::position(const LFHCALChipID& chip) const {
  return LFHCALCellPosition{
      10.0 * chip_x_cm(chip.moduleIDx, chip.module_type, chip.module_half),
      10.0 * chip_y_cm(chip.moduleIDy),
  };
}

int HcalEndcapPInsertCellIDDecoder::get(std::uint64_t cell_id, std::string_view field) const {
  if (field == "system") return decode_bits(cell_id, kInsertSystemOffset, kInsertSystemWidth);
  if (field == "side") return decode_bits(cell_id, kInsertSideOffset, kInsertSideWidth);
  if (field == "layer") return decode_bits(cell_id, kInsertLayerOffset, kInsertLayerWidth);
  if (field == "slice") return decode_bits(cell_id, kInsertSliceOffset, kInsertSliceWidth);
  if (field == "x") return decode_signed_bits(cell_id, kInsertXOffset, kInsertXWidth);
  if (field == "y") return decode_signed_bits(cell_id, kInsertYOffset, kInsertYWidth);
  throw std::invalid_argument("Unknown HcalEndcapPInsert cellID field");
}

HcalEndcapPInsertCellID HcalEndcapPInsertCellIDDecoder::cell(std::uint64_t cell_id) const {
  return HcalEndcapPInsertCellID{
      get(cell_id, "side"),
      get(cell_id, "layer"),
      get(cell_id, "slice"),
      get(cell_id, "x"),
      get(cell_id, "y"),
  };
}

}  // namespace rates
