#pragma once

#include "occupancy/include/shared.h"

#include <TDirectory.h>

#include <array>
#include <cstdint>
#include <unordered_map>

namespace rates::lfhcal_analysis::occupancy::central {

using ChannelPassMap = std::unordered_map<rates::LFHCALChannelID, std::uint64_t, rates::LFHCALChannelIDHash>;

struct Outputs {
  // Inclusive and origin-separated channel passes inside the central radius cut.
  std::array<ChannelPassMap, kMIPCoefficients.size()> all_passes;
  std::array<std::array<ChannelPassMap, kMIPCoefficients.size()>, 2> origin_passes;
};

void accumulate_event(Outputs& outputs, const EventChannels& event_channels);
void write_output(TDirectory* parent, const Outputs& outputs, std::uint64_t n_events);

}  // namespace rates::lfhcal_analysis::occupancy::central
