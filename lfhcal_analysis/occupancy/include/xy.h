#pragma once

#include "occupancy/include/shared.h"

#include <array>

namespace rates::lfhcal_analysis::occupancy::xy {

constexpr std::size_t kAllOrigins = rates::lfhcal_analysis::kNOriginFamilies;
constexpr std::size_t kNOriginGroups = kAllOrigins + 1;

struct Outputs {
  // Nominal-threshold channel maps for every origin and the inclusive group.
  std::array<std::array<LayerSum, kNReadoutLayers + 1>, kNOriginGroups> groups;
};

void accumulate_event(std::array<LayerSum, kNReadoutLayers + 1>& layer_sums,
                      std::array<TH1D*, kNChannelTypes>& channel_type_event_hists,
                      Outputs& outputs,
                      const EventChannels& passing_channels);
void write_output(TDirectory* parent, const Outputs& outputs, std::uint64_t n_events);

}  // namespace rates::lfhcal_analysis::occupancy::xy
