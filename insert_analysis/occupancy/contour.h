#pragma once

#include <TDirectory.h>

#include <vector>

namespace rates::insert_analysis::contour {

void write_output(TDirectory* parent,
                  const std::vector<double>& chip_rates_hz,
                  const std::vector<double>& mean_fired_channels,
                  double x_max,
                  double y_min,
                  double y_max,
                  const std::vector<double>& contour_rates_gbps,
                  const char* title_suffix);

}  // namespace rates::insert_analysis::contour
