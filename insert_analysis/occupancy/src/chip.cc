#include "chip.h"

#include "shared.h"

#include <TH1D.h>

#include <string>
#include <unordered_map>

namespace rates::insert_occupancy::chip {
namespace {
// ----------------------------------------------------------------------------------
// Constants
// ----------------------------------------------------------------------------------

// Transverse virtual chip size = 10 by 10 cm so it can cover four channels of 5 by 5 in a square.
constexpr double kChipSizeMM = 100.0;

// ----------------------------------------------------------------------------------
// Plot helper(s)
// ----------------------------------------------------------------------------------
// Called when writing virtual-chip rate histograms.
void write_chip_directory(TDirectory* parent,
                          const ChipSum& chip_sum,
                          std::uint64_t n_events) {
  const int segment = rates::insert_occupancy::kAllSegmentsIndex;
  const std::string dir_name = segment_dir_name(segment);
  auto* dir = parent->mkdir(dir_name.c_str());
  dir->cd();

  auto* hist_dir = dir->mkdir("hists");
  hist_dir->cd();
  chip_sum.h_hits_evt->Write();

  auto* h_hit_rate = new TH1D(
      "h_hit_rate",
      (segment_title(segment) + ";rate [Hz/virtual chip];virtual chips").c_str(),
      200,
      0.0,
      2.0e5);

  auto* h_data_rate = new TH1D(
      "h_data_rate",
      (segment_title(segment) + ";tail data rate [Gb/s/virtual chip];virtual chips").c_str(),
      200,
      0.0,
      0.05);

  // Convert each accumulated virtual chip into one full-sample rate entry.
  for (const auto& [chip, stats] : chip_sum.chips) {
    (void)chip;
    const double hit_rate = static_cast<double>(stats.total_hits) /
                            (static_cast<double>(n_events) * rates::insert_occupancy::kEventWindowSec);
    const double data_rate = stats.total_payload_bits /
                             (static_cast<double>(n_events) * rates::insert_occupancy::kEventWindowSec * 1.0e9);
    h_hit_rate->Fill(hit_rate);
    h_data_rate->Fill(data_rate);
  }

  h_hit_rate->Write();
  h_data_rate->Write();

  rates::insert_occupancy::draw_and_write(dir, h_hit_rate, "c_hit_rate", false, true);
  rates::insert_occupancy::draw_and_write(dir, h_data_rate, "c_data_rate", false, true);
  rates::insert_occupancy::draw_and_write(dir, chip_sum.h_hits_evt, "c_hits_evt", false, true);
}

}  // namespace

// Called once before the file loop to allocate chip-level products.
void init_chip_sum(ChipSum& chip_sum) {
  chip_sum.h_hits_evt = new TH1D(
      "h_hits_evt",
      (segment_title(rates::insert_occupancy::kAllSegmentsIndex) + ";hits/event/chip;Events").c_str(),
      200,
      0,
      2000);
}

// Called per event after the hit loop to accumulate chip-level occupancy.
void accumulate_event(ChipSum& chip_sum,
                      const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& channel_hit_counts,
                      const rates::InsertToLFHCALMapper& mapper) {
  std::unordered_map<rates::VirtualLFHCALChipID, int, rates::VirtualLFHCALChipIDHash> chip_active_channel_counts;
  int chip_hits_per_event = 0;

  // First collapse the event's fired virtual channels into per-chip active-channel counts.
  for (const auto& [channel, count] : channel_hit_counts) {
    (void)count;
    const auto chip = mapper.chip(channel);
    ++chip_active_channel_counts[chip]; // Used to track payload.
  }

  // Then convert each fired chip into occupancy and data-volume counters.
  for (const auto& [chip, active_channel_count] : chip_active_channel_counts) {
    const double payload_bits =
        (rates::insert_occupancy::kOverheadBits +
         rates::insert_occupancy::kBitsPerHit * static_cast<double>(active_channel_count)) *
        rates::insert_occupancy::kSamplesPerEvent;

    auto& stats = chip_sum.chips[chip];
    stats.x_mm = (static_cast<double>(chip.ix) + 0.5) * kChipSizeMM;
    stats.y_mm = (static_cast<double>(chip.iy) + 0.5) * kChipSizeMM;
    ++stats.total_hits;         // Count one fired event for this virtual chip.
    stats.total_payload_bits += payload_bits;  // Accumulate the chip payload from this event.

    ++chip_hits_per_event;
  }

  // Fill the per-event chip multiplicity histograms once the event is fully summed.
  chip_sum.h_hits_evt->Fill(chip_hits_per_event);
}

// Called once after the event loop to write virtual-chip products.
void write_output(TFile& output,
                  const ChipSum& chip_sum,
                  std::uint64_t n_events) {
  auto* parent = output.mkdir("chip");
  write_chip_directory(parent, chip_sum, n_events);
}

}  // namespace rates::insert_occupancy::chip
