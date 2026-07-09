#pragma once

#include "insert_to_lfhcal.h"

#include <TFile.h>
#include <TH1D.h>

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace rates::insert_occupancy::chip {

// ----------------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------------
struct ChipStats {
  double x_mm = 0.0;
  double y_mm = 0.0;
  std::uint64_t total_hits = 0;
  double total_bits = 0.0;
};

struct LayerAccum {
  std::unordered_map<rates::VirtualLFHCALChipID, ChipStats, rates::VirtualLFHCALChipIDHash> chips;
  TH1D* h_hits_evt = nullptr;
};

// ----------------------------------------------------------------------------------
// Interface
// ----------------------------------------------------------------------------------
void init_layer_accumulations(std::vector<LayerAccum>& layers);
void accumulate_event(std::vector<LayerAccum>& layers,
                      const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper);
void write_output(TFile& output,
                  const std::vector<LayerAccum>& layers,
                  std::uint64_t n_events);

}  // namespace rates::insert_occupancy::chip
