#pragma once

#include "occupancy/include/shared.h"

#include <TDirectory.h>

#include <array>
#include <cstdint>
#include <unordered_map>

namespace rates::lfhcal_analysis::occupancy::mip {

struct ThresholdSum {
  // Channel, chip, and payload totals accumulated at one MIP threshold.
  std::unordered_map<rates::LFHCALChannelID, std::uint64_t, rates::LFHCALChannelIDHash> channel_passes;
  std::unordered_map<rates::LFHCALChipID, std::uint64_t, rates::LFHCALChipIDHash> chip_passes;
  std::unordered_map<rates::LFHCALChipID, double, rates::LFHCALChipIDHash> chip_payload_bits;
  std::uint64_t total_active_channels = 0;
  std::uint64_t total_chip_instances = 0;
};

void accumulate_event(std::array<ThresholdSum, kMIPCoefficients.size()>& threshold_sums,
                      const EventChannels& event_channels);
void write_output(TDirectory* parent,
                  const std::array<ThresholdSum, kMIPCoefficients.size()>& threshold_sums,
                  std::uint64_t n_events);

}  // namespace rates::lfhcal_analysis::occupancy::mip
