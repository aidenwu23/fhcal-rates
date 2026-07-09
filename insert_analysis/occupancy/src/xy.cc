#include "xy.h"

#include <TH2D.h>

#include <string>

namespace rates::insert_occupancy::xy {
namespace {
// ----------------------------------------------------------------------------------
// Plot helper(s)
// ----------------------------------------------------------------------------------

// Called when writing one segment directory into the output file.
void write_segment_directory(TDirectory* parent,
                             const SegmentSum& segment_sum,
                             int segment,
                             std::uint64_t n_events) {
  const std::string dir_name = segment_dir_name(segment);
  auto* dir = parent->mkdir(dir_name.c_str());
  dir->cd();

  auto* hist_dir = dir->mkdir("hists");
  hist_dir->cd();
  segment_sum.h_hits_evt->Write();

  // Leave early when this segment never received any virtual channels.
  if (segment_sum.channels.empty()) {
    return;
  }

  // Build x and y binning from the channels that were actually filled.
  const auto axis_edges = make_segment_axis_edges(segment_sum);
  const auto& x_edges = axis_edges.x_edges;
  const auto& y_edges = axis_edges.y_edges;

  auto* h_avg = new TH2D(
      "h_avg",
      (segment_title(segment) + ";x [mm];y [mm];avg hits/event/virtual channel").c_str(),
      static_cast<int>(x_edges.size()) - 1,
      x_edges.data(),
      static_cast<int>(y_edges.size()) - 1,
      y_edges.data());

  auto* h_rate = new TH2D(
      "h_rate",
      (segment_title(segment) + ";x [mm];y [mm];rate [Hz/virtual channel]").c_str(),
      static_cast<int>(x_edges.size()) - 1,
      x_edges.data(),
      static_cast<int>(y_edges.size()) - 1,
      y_edges.data());

  // Fill one x-y bin per virtual channel using the accumulated full-sample totals.
  for (const auto& [channel, stats] : segment_sum.channels) {
    (void)channel;
    const double avg = static_cast<double>(stats.total_hits) / static_cast<double>(n_events);
    const double rate = static_cast<double>(stats.total_hits) / (static_cast<double>(n_events) * kEventWindowSec);
    h_avg->Fill(stats.x_mm, stats.y_mm, avg);
    h_rate->Fill(stats.x_mm, stats.y_mm, rate);
  }

  hist_dir->cd();
  h_avg->Write();
  h_rate->Write();

  draw_and_write(dir, h_avg, "c_avg", true, false);
  draw_and_write(dir, h_rate, "c_rate", true, false);
  draw_and_write(dir, segment_sum.h_hits_evt, "c_hits_evt", false, true);
}

}  // namespace

// Called per event after the hit loop to accumulate x-y occupancy.
void accumulate_event(std::vector<SegmentSum>& segment_sums,
                      const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& channel_hit_counts) {
  // Each event contributes to its own segment and to the inclusive summed-segment map.
  for (const auto& [channel, count] : channel_hit_counts) {
    const int segment = segment_index(channel.layer);
    if (segment < 0) continue;

    auto& segment_sum = segment_sums[segment];
    auto& stats = segment_sum.channels[channel];
    stats.x_mm = (static_cast<double>(channel.ix) + 0.5) * kVirtualCellSizeMM;
    stats.y_mm = (static_cast<double>(channel.iy) + 0.5) * kVirtualCellSizeMM;
    stats.total_hits += count;

    auto& all_segment_stats = segment_sums[kAllSegmentsIndex].channels[channel];
    all_segment_stats.x_mm = stats.x_mm;
    all_segment_stats.y_mm = stats.y_mm;
    all_segment_stats.total_hits += count;
  }
}

// Called once after the event loop to write x-y occupancy products.
void write_output(TFile& output,
                  const std::vector<SegmentSum>& segment_sums,
                  std::uint64_t n_events) {
  auto* parent = output.mkdir("xy");
  for (int segment = 0; segment <= kNSegments; ++segment) {
    write_segment_directory(parent, segment_sums[segment], segment, n_events);
  }
}

}  // namespace rates::insert_occupancy::xy
