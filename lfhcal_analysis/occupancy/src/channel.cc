#include "occupancy/include/channel.h"

#include <TH1D.h>

#include <string>

namespace rates::lfhcal_analysis::occupancy::channel {

void init(std::array<LayerSum, kNReadoutLayers + 1>& layer_sums) {
  // Create one fired-channel count histogram for each layer group.
  for (int layer = 0; layer <= kAllLayers; ++layer) {
    layer_sums[layer].h_channels_event = new TH1D(
        ("h_" + layer_name(layer) + "_channels_event").c_str(),
        (layer_title(layer) + ";fired channels/event;events").c_str(),
        200,
        0.0,
        2000.0);
  }
}

void write_output(TDirectory* parent,
                  const std::array<LayerSum, kNReadoutLayers + 1>& layer_sums,
                  std::uint64_t n_events) {
  // Create the channel hit-rate and data-rate output groups.
  auto* hit_dir = parent->mkdir("channel_hit_rate");
  auto* data_dir = parent->mkdir("channel_data_rate");
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;

  // Build rate distributions for each layer and the inclusive layer sum.
  for (int layer = 0; layer <= kAllLayers; ++layer) {
    const std::string name = layer_name(layer);

    hit_dir->cd();
    TH1D hit_rate(("h_" + name + "_hit_rate").c_str(),
                  (layer_title(layer) + ";rate [Hz/channel];channels").c_str(),
                  240,
                  0.0,
                  2.0e5);
    data_dir->cd();
    TH1D data_rate(("h_" + name + "_data_rate").c_str(),
                   (layer_title(layer) + ";data rate [Gb/s/channel];channels").c_str(),
                   240,
                   0.0,
                   0.03);

    // Convert each channel's passing-event count into hit and payload rates.
    for (const auto& [channel_id, stats] : layer_sums[layer].channels) {
      (void)channel_id;
      hit_rate.Fill(static_cast<double>(stats.passing_events) / total_time_sec);
      data_rate.Fill(static_cast<double>(stats.passing_events) * kBitsPerHit * kSamplesPerEvent /
                     (total_time_sec * 1.0e9));
    }

    // Write this layer's rate and event-multiplicity histograms.
    hit_dir->cd();
    hit_rate.Write();
    layer_sums[layer].h_channels_event->Write();
    data_dir->cd();
    data_rate.Write();
  }

  // Draw the inclusive event-multiplicity summary.
  draw_and_write(hit_dir, layer_sums[kAllLayers].h_channels_event, "c_channels_event", false, true);
}

}  // namespace rates::lfhcal_analysis::occupancy::channel
