#pragma once

#include "original_insert/include/shared.h"

#include <unordered_map>
#include <vector>

namespace rates::insert_analysis::original_insert::xy {

void accumulate_event(std::vector<SegmentSum>& segment_sums,
                      const OccupancyMode& mode,
                      const std::unordered_map<OriginalChannelID, int, OriginalChannelIDHash>& channel_hit_counts);

}  // namespace rates::insert_analysis::original_insert::xy
