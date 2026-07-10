#pragma once

namespace rates::insert_analysis {

struct EventHit {
  // Decoded insert cell coordinates and energy used by both occupancy layouts.
  int layer = 0;
  int side = 0;
  double x_mm = 0.0;
  double y_mm = 0.0;
  double energy_gev = 0.0;
};

}  // namespace rates::insert_analysis
