#pragma once

#include <TFile.h>

#include "event_hit.h"
#include "insert_to_lfhcal.h"
#include "original_insert/include/chip.h"
#include "original_insert/include/mip.h"
#include "original_insert/include/shared.h"
#include "original_insert/include/side.h"

#include <array>
#include <vector>

namespace rates::insert_analysis::original_insert {

struct Outputs {
  // Keep all products for one inner-ring treatment together.
  OccupancyMode mode;
  std::vector<SegmentSum> segment_sums;
  chip::ChipSum chip_sum;
  std::array<mip::ChannelThresholdSum, 16> channel_threshold_sums;
  std::array<mip::ChannelDataThresholdSum, 16> channel_data_threshold_sums;
  std::array<mip::ChipThresholdSum, 16> chip_threshold_sums;
  std::array<mip::ChipDataThresholdSum, 16> chip_data_threshold_sums;
  std::array<std::array<side::ThresholdSum, 16>, 2> side_threshold_sums;
};

void init_outputs(Outputs& outputs);
void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper);
void write_output(TFile& output,
                  const Outputs& outputs,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::original_insert
