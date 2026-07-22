#include "occupancy/include/multiplicity.h"

#include <string>

namespace rates::lfhcal_analysis::occupancy::multiplicity {

void init(Outputs& outputs) {
  // Create one hit-multiplicity histogram per readout layer.
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    outputs.layer_hists[layer] = new TH1D(
        ("h_channel_multiplicity_layer" + std::to_string(layer)).c_str(),
        ("LFHCAL channel multiplicity per event, layer " + std::to_string(layer + 1) +
         ";hits/channel/event;channels").c_str(),
        100,
        0.0,
        100.0);
  }
}

void accumulate_event(Outputs& outputs, const EventChannels& passing_channels) {
  // Fill the number of simulated hits summed into each passing channel.
  for (const auto& event_channel : passing_channels) {
    outputs.layer_hists[event_channel.channel.rlayerz]->Fill(event_channel.hit_count);
  }
}

void write_output(TDirectory* parent, const Outputs& outputs) {
  // Write every layer's multiplicity histogram and canvas.
  auto* directory = parent->mkdir("channel_multiplicity");
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    directory->cd();
    outputs.layer_hists[layer]->Write();
    draw_and_write(directory,
                   outputs.layer_hists[layer],
                   ("c_layer" + std::to_string(layer)).c_str(),
                   false,
                   true);
  }
}

}  // namespace rates::lfhcal_analysis::occupancy::multiplicity
