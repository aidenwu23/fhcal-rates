#pragma once

#include "insert_to_lfhcal.h"
#include "lfhcal_tiles/include/shared.h"

#include <TDirectory.h>
#include <TH1D.h>

#include <array>
#include <cstdint>
#include <unordered_map>

namespace rates::insert_analysis::lfhcal_tiles::chip {

struct ChipStats {
  std::uint64_t total_hits = 0;
  double total_payload_bits = 0.0;
};

struct ChipSum {
  std::unordered_map<rates::VirtualLFHCALPizzaChipID,
                     ChipStats,
                     rates::VirtualLFHCALPizzaChipIDHash>
      chips;
  TH1D* h_hits_evt = nullptr;
};

using ChannelHitCountMap =
    std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>;

void init_chip_sum(ChipSum& chip_sum, const OccupancyMode& mode);
void accumulate_event(ChipSum& chip_sum,
                      const std::array<ChannelHitCountMap, 2>& side_channel_hit_counts,
                      const rates::InsertToLFHCALMapper& mapper);
void write_output(TDirectory* parent,
                  const OccupancyMode& mode,
                  const ChipSum& chip_sum,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::lfhcal_tiles::chip
