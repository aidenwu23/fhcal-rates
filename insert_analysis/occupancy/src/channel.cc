#include "channel.h"

#include <TH1D.h>

#include <string>

namespace rates::insert_occupancy::channel {
namespace {

// ----------------------------------------------------------------------------------
// Plot helper(s)
// ----------------------------------------------------------------------------------

// Called when writing one layer of virtual-channel rate histograms.
void write_layer_directory(TDirectory* parent,
                           const LayerAccum& layer_accum,
                           int layer,
                           std::uint64_t n_events) {
  const std::string dir_name = segment_dir_name(layer);
  auto* dir = parent->mkdir(dir_name.c_str());
  dir->cd();

  auto* hist_dir = dir->mkdir("hists");
  hist_dir->cd();

  auto* h_hit_rate = new TH1D(
      "h_hit_rate",
      (segment_title(layer) + ";rate [Hz/virtual channel];virtual channels").c_str(),
      200,
      0.0,
      5.0e4);

  auto* h_data_rate = new TH1D(
      "h_data_rate",
      (segment_title(layer) + ";data rate [Gb/s/virtual channel];virtual channels").c_str(),
      200,
      0.0,
      0.01);

  // Convert each accumulated virtual channel into one full-sample rate entry.
  for (const auto& [channel, stats] : layer_accum.channels) {
    (void)channel;
    const double hit_rate = static_cast<double>(stats.total_hits) / (static_cast<double>(n_events) * kEventWindowSec);
    const double data_rate =
        (static_cast<double>(stats.total_hits) * kBitsPerHit * kSamplesPerEvent) /
        (static_cast<double>(n_events) * kEventWindowSec * 1.0e9);
    h_hit_rate->Fill(hit_rate);
    h_data_rate->Fill(data_rate);
  }

  h_hit_rate->Write();
  h_data_rate->Write();

  draw_and_write(dir, h_hit_rate, "c_hit_rate", false, true);
  draw_and_write(dir, h_data_rate, "c_data_rate", false, true);
}

}  // namespace

// Called once after the event loop to write virtual-channel products.
void write_output(TFile& output,
                  const std::vector<LayerAccum>& layers,
                  std::uint64_t n_events) {
  auto* parent = output.mkdir("channel");

  // Write one directory per segment plus one inclusive summed-segment directory.
  for (int layer = 0; layer <= kAllSegmentsIndex; ++layer) {
    write_layer_directory(parent, layers[layer], layer, n_events);
  }
}

}  // namespace rates::insert_occupancy::channel
