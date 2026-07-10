#pragma once

#include <TFile.h>

#include "event_hit.h"
#include "insert_to_lfhcal.h"
#include "physical_layers/include/mip.h"
#include "physical_layers/include/shared.h"

#include <array>
#include <vector>

namespace rates::insert_analysis::physical_layers {

struct Outputs {
  // Keep all products for one inner-ring treatment together.
  OccupancyMode mode;
  std::vector<SegmentSum> segment_sums;
  std::array<mip::ChannelThresholdSum, 16> channel_threshold_sums;
  std::array<mip::ChannelDataThresholdSum, 16> channel_data_threshold_sums;
};

void init_outputs(Outputs& outputs);
void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper);
void write_output(TFile& output,
                  const Outputs& outputs,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::physical_layers
