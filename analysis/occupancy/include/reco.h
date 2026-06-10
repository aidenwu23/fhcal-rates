#pragma once

#include <TH1D.h>

#include <podio/Frame.h>

#include "decode_cell_id.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace br::occupancy {

constexpr int kNLayers = 7;
constexpr int kAllLayersIndex = kNLayers;

// Stats per channel (cells sharing transverse location in a readout layer).
struct ChannelStats {
  std::unordered_set<std::uint64_t> raw_cell_ids; // Keep track of which cells belong in this channel.

  double x_mm = 0.0;
  double y_mm = 0.0;
  bool has_position = false;

  std::uint64_t total_hits = 0;

  int max_hits_event = 0; // Keep track of which event had the most hits.

  void set_position(double x, double y);
  double x() const;
  double y() const;
  double r() const; // radial distance from center.
};

// Hold one layer's data.
struct LayerAccum {

  // For each unique channelID, it stores the statistics and also the ID hash.
  std::unordered_map<br::LFHCALChannelID, ChannelStats, br::LFHCALChannelIDHash> channels;
  TH1D* h_hits_evt = nullptr;
};

void init_reco_layers(std::vector<LayerAccum>& layers);
bool process_reco_event(const podio::Frame& frame,
                        const br::LFHCALCellIDDecoder& decoder,
                        double threshold_geV,
                        std::vector<LayerAccum>& layers);

}  // namespace br::occupancy
