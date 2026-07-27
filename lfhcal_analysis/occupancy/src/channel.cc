#include "occupancy/include/channel.h"

#include <TH1D.h>

#include <string>

namespace rates::lfhcal_analysis::occupancy::channel {

void init(std::array<TH1D*, kNChannelTypes>& channel_type_event_hists) {
  // Create one event-wise fired-channel histogram for each channel type.
  for (int type = 0; type < kNChannelTypes; ++type) {
    channel_type_event_hists[type] = new TH1D(
        ("h_" + channel_type_name(type) + "_channels_event").c_str(),
        (channel_type_title(type) + ";fired channels/event;events").c_str(),
        200,
        0.0,
        2000.0);
  }
}

void write_output(TDirectory* parent,
                  const std::array<LayerSum, kNReadoutLayers + 1>& layer_sums,
                  const std::array<TH1D*, kNChannelTypes>& channel_type_event_hists,
                  std::uint64_t n_events) {
  // Create the channel hit-rate and data-rate output groups.
  auto* hit_dir = parent->mkdir("channel_hit_rate");
  auto* data_dir = parent->mkdir("channel_data_rate");
  auto* hit_hist_dir = hit_dir->mkdir("hists");
  auto* data_hist_dir = data_dir->mkdir("hists");
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;

  // Build pooled rate distributions for 5-tile and 10-tile channels.
  for (int type = 0; type < kNChannelTypes; ++type) {
    const std::string name = channel_type_name(type);
    const int first_layer = type == 0 ? 0 : 2;
    const int last_layer = type == 0 ? 1 : 6;

    hit_hist_dir->cd();
    TH1D hit_rate(("h_" + name + "_hit_rate").c_str(),
                  (channel_type_title(type) + ";rate [Hz/channel];channels").c_str(),
                  240,
                  0.0,
                  2.0e5);
    data_hist_dir->cd();
    TH1D data_rate(("h_" + name + "_data_rate").c_str(),
                   (channel_type_title(type) + ";data rate [Gb/s/channel];channels").c_str(),
                   240,
                   0.0,
                   0.03);

    // Pool channels from every layer with the same longitudinal tile count.
    for (int layer = first_layer; layer <= last_layer; ++layer) {
      for (const auto& [channel_id, stats] : layer_sums[layer].channels) {
        (void)channel_id;
        hit_rate.Fill(static_cast<double>(stats.passing_events) / total_time_sec);
        data_rate.Fill(static_cast<double>(stats.passing_events) * kBitsPerHit * kSamplesPerEvent /
                       (total_time_sec * 1.0e9));
      }
    }

    // Write this channel type's histograms and canvases.
    hit_hist_dir->cd();
    hit_rate.Write();
    channel_type_event_hists[type]->Write();
    draw_and_write(hit_dir,
                   &hit_rate,
                   ("c_" + name + "_hit_rate").c_str(),
                   false,
                   true);
    draw_and_write(hit_dir,
                   channel_type_event_hists[type],
                   ("c_" + name + "_channels_event").c_str(),
                   false,
                   true);
    data_hist_dir->cd();
    data_rate.Write();
    draw_and_write(data_dir,
                   &data_rate,
                   ("c_" + name + "_data_rate").c_str(),
                   false,
                   true);
  }

}

}  // namespace rates::lfhcal_analysis::occupancy::channel
