#pragma once

#include "occupancy/include/channel.h"
#include "occupancy/include/central.h"
#include "occupancy/include/chip.h"
#include "occupancy/include/mip.h"
#include "occupancy/include/multiplicity.h"
#include "occupancy/include/radius.h"
#include "occupancy/include/xy.h"

#include <TFile.h>

#include <array>
#include <cstdint>

namespace rates::lfhcal_analysis::occupancy {

struct Outputs {
  // Full-sample accumulators for every occupancy output family.
  std::array<LayerSum, kNReadoutLayers + 1> layer_sums;
  std::array<TH1D*, kNChannelTypes> channel_type_event_hists{};
  xy::Outputs spatial_outputs;
  chip::ChipSum chip_sum;
  central::Outputs central_outputs;
  std::array<mip::ThresholdSum, kMIPCoefficients.size()> mip_threshold_sums;
  multiplicity::Outputs multiplicity_outputs;
  std::array<radius::ThresholdSum, radius::kThresholds.size()> radius_threshold_sums;
};

void init_outputs(Outputs& outputs);
void accumulate_event(Outputs& outputs,
                      const EventChannels& event_channels,
                      const rates::LFHCALCellIDDecoder& decoder);
void write_output(TFile& output,
                  const Outputs& outputs,
                  const rates::LFHCALCellIDDecoder& decoder,
                  std::uint64_t n_events);

}  // namespace rates::lfhcal_analysis::occupancy
