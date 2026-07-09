#pragma once

#include "shared.h"

#include <TFile.h>

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace rates::insert_occupancy::xy {

void accumulate_event(std::vector<SegmentSum>& segment_sums,
                      const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& channel_hit_counts);

void write_output(TFile& output,
                  const std::vector<SegmentSum>& segment_sums,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::xy
