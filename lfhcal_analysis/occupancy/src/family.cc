#include "occupancy/include/family.h"

#include "occupancy/include/xy.h"

namespace rates::lfhcal_analysis::occupancy {

void init_outputs(Outputs& outputs) {
  // Initialize every histogram-based output family.
  channel::init(outputs.channel_type_event_hists);
  chip::init(outputs.chip_sum);
  multiplicity::init(outputs.multiplicity_outputs);
}

void accumulate_event(Outputs& outputs,
                      const EventChannels& event_channels,
                      const rates::LFHCALCellIDDecoder& decoder) {
  // Apply the nominal channel threshold.
  EventChannels passing_channels;
  passing_channels.reserve(event_channels.size());

  // Select channels that pass the nominal threshold.
  for (const auto& event_channel : event_channels) {
    const int layer = event_channel.channel.rlayerz;
    if (layer < 0 || layer >= kNReadoutLayers) continue;
    if (event_channel.energy_gev <= kNominalThresholdMIP * rates::mip_energy_gev(layer)) continue;
    passing_channels.push_back(event_channel);
  }

  // Fill fixed-threshold channel, chip, and multiplicity products.
  xy::accumulate_event(outputs.layer_sums,
                       outputs.channel_type_event_hists,
                       outputs.spatial_outputs,
                       passing_channels);
  chip::accumulate_event(outputs.chip_sum, passing_channels, decoder);
  multiplicity::accumulate_event(outputs.multiplicity_outputs, passing_channels);
  
  // Fill threshold-scan, radial, and central-region products.
  mip::accumulate_event(outputs.mip_threshold_sums, event_channels);
  radius::accumulate_event(outputs.radius_threshold_sums, event_channels);
  central::accumulate_event(outputs.central_outputs, event_channels);
}

void write_output(TFile& output,
                  const Outputs& outputs,
                  const rates::LFHCALCellIDDecoder& decoder,
                  std::uint64_t n_events) {
  // Write every analysis family beneath the occupancy directory.
  auto* occupancy_dir = output.mkdir("occupancy");
  channel::write_output(occupancy_dir, outputs.layer_sums, outputs.channel_type_event_hists, n_events);
  xy::write_output(occupancy_dir, outputs.spatial_outputs, n_events);
  chip::write_output(occupancy_dir, outputs.chip_sum, n_events);
  multiplicity::write_output(occupancy_dir, outputs.multiplicity_outputs);
  mip::write_output(occupancy_dir, outputs.mip_threshold_sums, n_events);
  radius::write_output(occupancy_dir, outputs.radius_threshold_sums, decoder, n_events);
  central::write_output(occupancy_dir, outputs.central_outputs, n_events);
}

}  // namespace rates::lfhcal_analysis::occupancy
