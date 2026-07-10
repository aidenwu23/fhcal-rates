#pragma once

#include "full_lfhcal/include/shared.h"

#include <unordered_map>
#include <vector>

namespace rates::insert_analysis::full_lfhcal::xy {

void accumulate_event(std::vector<SegmentSum>& segment_sums,
                      const OccupancyMode& mode,
                      const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& channel_hit_counts);

}  // namespace rates::insert_analysis::full_lfhcal::xy
