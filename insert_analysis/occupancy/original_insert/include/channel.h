#pragma once

#include "original_insert/include/shared.h"

#include <TDirectory.h>

#include <cstdint>
#include <vector>

namespace rates::insert_analysis::original_insert::channel {

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void write_output(TDirectory* parent,
                  const OccupancyMode& mode,
                  const std::vector<SegmentSum>& segment_sums,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::original_insert::channel
