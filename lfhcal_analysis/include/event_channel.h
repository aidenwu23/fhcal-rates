#pragma once

#include "decode_cell_id.h"

#include <array>

namespace rates::lfhcal_analysis {

struct EventChannel {
  rates::LFHCALChannelID channel;
  rates::LFHCALChipID chip;
  double x_mm = 0.0;
  double y_mm = 0.0;
  double energy_gev = 0.0;
  std::array<double, 2> energy_by_origin{};
  int hit_count = 0;
};

}  // namespace rates::lfhcal_analysis
