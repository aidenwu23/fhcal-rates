#include "lfhcal_segments/include/chip.h"

#include "lfhcal_segments/include/shared.h"

#include <TH1D.h>

#include <string>
#include <unordered_map>

namespace rates::insert_analysis::lfhcal_segments::chip {
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
void write_chip_hit_rate_directory(TDirectory* parent,
                          const OccupancyMode& mode,
                          const ChipSum& chip_sum,
                          std::uint64_t n_events) {
  const int group = rates::insert_analysis::lfhcal_segments::all_groups_index(mode);
  const std::string dir_name = group_dir_name(mode, group);
  auto* dir = parent->mkdir(dir_name.c_str());
  dir->cd();

  auto* hist_dir = dir->mkdir("hists");
  hist_dir->cd();
  chip_sum.h_hits_evt->Write();

  auto* h_hit_rate = new TH1D(
      "h_hit_rate",
      (group_title(mode, group) + ";rate [Hz/virtual chip];virtual chips").c_str(),
      200,
      0.0,
      2.0e5);

  for (const auto& [chip, stats] : chip_sum.chips) {
    (void)chip;
    const double hit_rate = static_cast<double>(stats.total_hits) /
                            (static_cast<double>(n_events) * rates::insert_analysis::lfhcal_segments::kEventWindowSec);
    h_hit_rate->Fill(hit_rate);
  }

  h_hit_rate->Write();
  rates::insert_analysis::lfhcal_segments::draw_and_write(dir, h_hit_rate, "c_hit_rate", false, true);
  rates::insert_analysis::lfhcal_segments::draw_and_write(dir, chip_sum.h_hits_evt, "c_hits_evt", false, true);
}

void write_chip_data_rate_directory(TDirectory* parent,
                                    const OccupancyMode& mode,
                                    const ChipSum& chip_sum,
                                    std::uint64_t n_events) {
  const int group = rates::insert_analysis::lfhcal_segments::all_groups_index(mode);
  const std::string dir_name = group_dir_name(mode, group);
  auto* dir = parent->mkdir(dir_name.c_str());
  dir->cd();

  auto* hist_dir = dir->mkdir("hists");
  hist_dir->cd();

  auto* h_data_rate = new TH1D(
      "h_data_rate",
      (group_title(mode, group) + ";tail data rate [Gb/s/virtual chip];virtual chips").c_str(),
      200,
      0.0,
      0.05);

  for (const auto& [chip, stats] : chip_sum.chips) {
    (void)chip;
    const double data_rate = stats.total_payload_bits /
                             (static_cast<double>(n_events) * rates::insert_analysis::lfhcal_segments::kEventWindowSec * 1.0e9);
    h_data_rate->Fill(data_rate);
  }

  h_data_rate->Write();
  rates::insert_analysis::lfhcal_segments::draw_and_write(dir, h_data_rate, "c_data_rate", false, true);
}

}  // namespace

// Called once before the file loop to allocate chip-level products.
void init_chip_sum(ChipSum& chip_sum, const OccupancyMode& mode) {
  // This histogram records the number of fired virtual chips in each event.
  chip_sum.h_hits_evt = new TH1D(
      "h_hits_evt",
      (group_title(mode, all_groups_index(mode)) + ";hits/event/chip;Events").c_str(),
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
    ++chip_active_channel_counts[chip]; // Count this fired channel in the chip payload.
  }

  // Then convert each fired chip into occupancy and data-volume counters.
  for (const auto& [chip, active_channel_count] : chip_active_channel_counts) {
    const double payload_bits =
        (rates::insert_analysis::lfhcal_segments::kOverheadBits +
         rates::insert_analysis::lfhcal_segments::kBitsPerHit * static_cast<double>(active_channel_count)) *
        rates::insert_analysis::lfhcal_segments::kSamplesPerEvent;

    auto& stats = chip_sum.chips[chip];
    stats.x_mm = (static_cast<double>(chip.ix) + 0.5) * kChipSizeMM;
    stats.y_mm = (static_cast<double>(chip.iy) + 0.5) * kChipSizeMM;
    ++stats.total_hits; // Count one fired event for this virtual chip.
    stats.total_payload_bits += payload_bits;  // Accumulate the chip payload from this event.

    ++chip_hits_per_event;
  }

  // Fill the per-event chip multiplicity histograms once the event is fully summed.
  chip_sum.h_hits_evt->Fill(chip_hits_per_event);
}

// Called once after the event loop to write virtual-chip products.
void write_output(TDirectory* parent,
                  const OccupancyMode& mode,
                  const ChipSum& chip_sum,
                  std::uint64_t n_events) {
  // Keep hit-rate and data-rate products in separate directory trees.
  auto* hit_dir = parent->mkdir("chip_hit_rate");
  auto* data_dir = parent->mkdir("chip_data_rate");
  write_chip_hit_rate_directory(hit_dir, mode, chip_sum, n_events);
  write_chip_data_rate_directory(data_dir, mode, chip_sum, n_events);
}

}  // namespace rates::insert_analysis::lfhcal_segments::chip
