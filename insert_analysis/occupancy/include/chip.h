#pragma once

#include "insert_to_lfhcal.h"

#include <TFile.h>
#include <TH1D.h>

#include <cstdint>
#include <unordered_map>

namespace rates::insert_occupancy::chip {

// ----------------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------------
struct ChipStats {
  double x_mm = 0.0;
  double y_mm = 0.0;
  std::uint64_t total_hits = 0;
  double total_payload_bits = 0.0;
};

struct ChipSum {
  std::unordered_map<rates::VirtualLFHCALChipID, ChipStats, rates::VirtualLFHCALChipIDHash> chips;
  TH1D* h_hits_evt = nullptr;
};

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void init_chip_sum(ChipSum& chip_sum);
void accumulate_event(ChipSum& chip_sum,
                      const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& channel_hit_counts,
                      const rates::InsertToLFHCALMapper& mapper);
void write_output(TFile& output,
                  const ChipSum& chip_sum,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::chip
