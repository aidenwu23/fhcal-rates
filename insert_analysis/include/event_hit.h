#pragma once

#include <cstdint>

namespace rates::insert_analysis {

struct EventHit {
  // Original insert cell identity, coordinates, and energy shared by all layouts.
  std::uint64_t cell_id = 0;
  int layer = 0;
  int side = 0;
  double x_mm = 0.0;
  double y_mm = 0.0;
  double energy_gev = 0.0;
};

}  // namespace rates::insert_analysis
