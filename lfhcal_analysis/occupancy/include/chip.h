#pragma once

#include "occupancy/include/shared.h"

#include <TDirectory.h>
#include <TH1D.h>

#include <cstdint>
#include <unordered_map>

namespace rates::lfhcal_analysis::occupancy::chip {

struct ChipStats {
  // Full-sample position, occupancy, multiplicity, and payload totals for one chip.
  double x_mm = 0.0;
  double y_mm = 0.0;
  std::uint64_t passing_events = 0;
  std::uint64_t total_active_channels = 0;
  double total_payload_bits = 0.0;
};

struct ChipSum {
  // Per-chip totals and event-level chip multiplicity histograms.
  std::unordered_map<rates::LFHCALChipID, ChipStats, rates::LFHCALChipIDHash> chips;
  TH1D* h_chips_event = nullptr;
  TH1D* h_active_channels = nullptr;
};

void init(ChipSum& chip_sum);
void accumulate_event(ChipSum& chip_sum,
                      const EventChannels& passing_channels,
                      const rates::LFHCALCellIDDecoder& decoder);
void write_output(TDirectory* parent, const ChipSum& chip_sum, std::uint64_t n_events);

}  // namespace rates::lfhcal_analysis::occupancy::chip
