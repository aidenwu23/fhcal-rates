#pragma once

#include <TDirectory.h>

#include <vector>

namespace rates::lfhcal_analysis::occupancy::contour {

void write_output(TDirectory* parent,
                  const std::vector<double>& chip_rates_hz,
                  const std::vector<double>& mean_fired_channels);

}  // namespace rates::lfhcal_analysis::occupancy::contour
