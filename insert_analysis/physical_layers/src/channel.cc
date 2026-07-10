#include "physical_layers/include/channel.h"

#include <TH1D.h>

#include <string>

namespace rates::insert_analysis::physical_layers::channel {
namespace {

// ----------------------------------------------------------------------------------
// Plot helper(s)
// ----------------------------------------------------------------------------------

void write_hit_rate_directory(TDirectory* parent,
                             const SegmentSum& segment_sum,
                             const OccupancyMode& mode,
                             int group,
                             std::uint64_t n_events) {
  // Keep every layer histogram directly in the channel-rate directory.
  const std::string group_name = group_dir_name(mode, group);
  const std::string hist_name = "h_" + group_name + "_hit_rate";
  parent->cd();

  auto* h_hit_rate = new TH1D(
      hist_name.c_str(),
      (group_title(mode, group) + ";rate [Hz/virtual channel];virtual channels").c_str(),
      200,
      0.0,
      5.0e4);

  // Convert each full-sample channel count into a rate using the event window.
  for (const auto& [channel, stats] : segment_sum.channels) {
    (void)channel;
    const double hit_rate = static_cast<double>(stats.total_hits) / (static_cast<double>(n_events) * kEventWindowSec);
    h_hit_rate->Fill(hit_rate);
  }

  h_hit_rate->Write();
}

void write_data_rate_directory(TDirectory* parent,
                               const SegmentSum& segment_sum,
                               const OccupancyMode& mode,
                               int group,
                               std::uint64_t n_events) {
  // Keep every layer histogram directly in the channel-rate directory.
  const std::string group_name = group_dir_name(mode, group);
  const std::string hist_name = "h_" + group_name + "_data_rate";
  parent->cd();

  auto* h_data_rate = new TH1D(
      hist_name.c_str(),
      (group_title(mode, group) + ";data rate [Gb/s/virtual channel];virtual channels").c_str(),
      200,
      0.0,
      0.01);

  // Convert hit multiplicity into the corresponding per-channel data rate.
  for (const auto& [channel, stats] : segment_sum.channels) {
    (void)channel;
    const double data_rate =
        (static_cast<double>(stats.total_hits) * kBitsPerHit * kSamplesPerEvent) /
        (static_cast<double>(n_events) * kEventWindowSec * 1.0e9);
    h_data_rate->Fill(data_rate);
  }

  h_data_rate->Write();
}

}  // namespace

// Called once after the event loop to write virtual-channel products.
void write_output(TDirectory* parent,
                  const OccupancyMode& mode,
                  const std::vector<SegmentSum>& segment_sums,
                  std::uint64_t n_events) {
  auto* hit_dir = parent->mkdir("channel_hit_rate");
  auto* data_dir = parent->mkdir("channel_data_rate");

  // Write every physical layer and the inclusive layer sum.
  for (int group = 0; group <= group_count(mode); ++group) {
    write_hit_rate_directory(hit_dir, segment_sums[group], mode, group, n_events);
    write_data_rate_directory(data_dir, segment_sums[group], mode, group, n_events);
  }
}

}  // namespace rates::insert_analysis::physical_layers::channel
