#include "occupancy/include/multiplicity.h"

#include <string>

namespace rates::lfhcal_analysis::occupancy::multiplicity {

void init(Outputs& outputs) {
  // Create one hit-multiplicity histogram per channel type.
  for (int type = 0; type < kNChannelTypes; ++type) {
    outputs.channel_type_hists[type] = new TH1D(
        ("h_channel_multiplicity_" + channel_type_name(type)).c_str(),
        (channel_type_title(type) + ";hits/channel/event;channels").c_str(),
        100,
        0.0,
        100.0);
  }
}

void accumulate_event(Outputs& outputs, const EventChannels& passing_channels) {
  // Fill the number of simulated hits summed into each passing channel.
  for (const auto& event_channel : passing_channels) {
    outputs.channel_type_hists[channel_type(event_channel.channel.rlayerz)]->Fill(event_channel.hit_count);
  }
}

void write_output(TDirectory* parent, const Outputs& outputs) {
  // Write every layer's multiplicity histogram and canvas.
  auto* directory = parent->mkdir("channel_multiplicity");
  auto* hist_directory = directory->mkdir("hists");
  for (int type = 0; type < kNChannelTypes; ++type) {
    hist_directory->cd();
    outputs.channel_type_hists[type]->Write();
    draw_and_write(directory,
                   outputs.channel_type_hists[type],
                   ("c_" + channel_type_name(type)).c_str(),
                   false,
                   true);
  }
}

}  // namespace rates::lfhcal_analysis::occupancy::multiplicity
