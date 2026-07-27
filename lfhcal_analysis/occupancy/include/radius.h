#pragma once

#include "occupancy/include/shared.h"

#include <TDirectory.h>

#include <array>
#include <cstdint>
#include <unordered_map>

namespace rates::lfhcal_analysis::occupancy::radius {

constexpr std::array<double, 3> kThresholds = {0.1, 0.3, 0.5};

struct ThresholdSum {
  std::unordered_map<rates::LFHCALChipID, double, rates::LFHCALChipIDHash> payload_bits;
};

void accumulate_event(std::array<ThresholdSum, kThresholds.size()>& threshold_sums,
                      const EventChannels& event_channels);
void write_output(TDirectory* parent,
                  const std::array<ThresholdSum, kThresholds.size()>& threshold_sums,
                  const rates::LFHCALCellIDDecoder& decoder,
                  std::uint64_t n_events);

}  // namespace rates::lfhcal_analysis::occupancy::radius
