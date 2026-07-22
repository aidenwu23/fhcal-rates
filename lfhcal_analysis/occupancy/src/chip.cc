#include "occupancy/include/chip.h"

#include "contour.h"

#include <TH1D.h>

#include <unordered_map>
#include <vector>

namespace rates::lfhcal_analysis::occupancy::chip {

void init(ChipSum& chip_sum) {
  // Create the event-level chip and channel multiplicity histograms.
  chip_sum.h_chips_event = new TH1D(
      "h_chips_event",
      "LFHCAL fired chips per event;fired chips/event;events",
      200,
      0.0,
      2000.0);
  chip_sum.h_active_channels = new TH1D(
      "h_active_channels",
      "LFHCAL fired channels per active chip-event;fired channels/chip/event;chip-events",
      29,
      -0.5,
      28.5);
}

void accumulate_event(ChipSum& chip_sum,
                      const EventChannels& passing_channels,
                      const rates::LFHCALCellIDDecoder& decoder) {
  // Group this event's passing channels by chip.
  std::unordered_map<rates::LFHCALChipID, int, rates::LFHCALChipIDHash> active_channels;
  for (const auto& event_channel : passing_channels) ++active_channels[event_channel.chip];

  // Accumulate active events, channel counts, and payload bits for each fired chip.
  for (const auto& [chip_id, active_count] : active_channels) {
    chip_sum.h_active_channels->Fill(active_count);
    auto& stats = chip_sum.chips[chip_id];
    const auto position = decoder.position(chip_id);
    stats.x_mm = position.x_mm;
    stats.y_mm = position.y_mm;
    ++stats.passing_events;
    stats.total_active_channels += static_cast<std::uint64_t>(active_count);
    stats.total_payload_bits +=
        (kOverheadBits + kBitsPerHit * static_cast<double>(active_count)) * kSamplesPerEvent;
  }
  // Record how many chips fired in this event.
  chip_sum.h_chips_event->Fill(active_channels.size());
}

void write_output(TDirectory* parent, const ChipSum& chip_sum, std::uint64_t n_events) {
  // Create chip hit-rate and data-rate distributions.
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  auto* hit_dir = parent->mkdir("chip_hit_rate");
  auto* data_dir = parent->mkdir("chip_data_rate");

  TH1D hit_rate("h_hit_rate", "LFHCAL chip hit rate;rate [Hz/chip];chips", 240, 0.0, 2.0e5);
  TH1D data_rate("h_data_rate", "LFHCAL chip data rate;data rate [Gb/s/chip];chips", 240, 0.0, 0.1);
  std::vector<double> chip_rates_hz;
  std::vector<double> mean_fired_channels;

  // Convert each chip's totals into rates and conditional channel multiplicity.
  for (const auto& [chip_id, stats] : chip_sum.chips) {
    (void)chip_id;
    const double rate_hz = static_cast<double>(stats.passing_events) / total_time_sec;
    hit_rate.Fill(rate_hz);
    data_rate.Fill(stats.total_payload_bits / (total_time_sec * 1.0e9));
    chip_rates_hz.push_back(rate_hz);
    mean_fired_channels.push_back(
        static_cast<double>(stats.total_active_channels) / static_cast<double>(stats.passing_events));
  }

  // Write the chip histograms and canvases.
  hit_dir->cd();
  hit_rate.Write();
  chip_sum.h_chips_event->Write();
  chip_sum.h_active_channels->Write();
  draw_and_write(hit_dir, &hit_rate, "c_hit_rate", false, true);
  data_dir->cd();
  data_rate.Write();
  draw_and_write(data_dir, &data_rate, "c_data_rate", false, true);

  // Write fired-channel multiplicity versus chip rate.
  auto* contour_dir = parent->mkdir("chip_fired_channels_vs_rate");
  contour::write_output(contour_dir, chip_rates_hz, mean_fired_channels);
}

}  // namespace rates::lfhcal_analysis::occupancy::chip
