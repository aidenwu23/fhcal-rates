#pragma once

#include "occupancy/include/shared.h"

#include <TDirectory.h>

#include <array>

namespace rates::lfhcal_analysis::occupancy::multiplicity {

struct Outputs {
  std::array<TH1D*, kNReadoutLayers> layer_hists{};
};

void init(Outputs& outputs);
void accumulate_event(Outputs& outputs, const EventChannels& passing_channels);
void write_output(TDirectory* parent, const Outputs& outputs);

}  // namespace rates::lfhcal_analysis::occupancy::multiplicity
