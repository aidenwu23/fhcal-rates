#include "lfhcal_tiles/include/chip.h"

#include "lfhcal_tiles/include/shared.h"

#include <TH1D.h>

#include <string>
#include <unordered_map>

namespace rates::insert_analysis::lfhcal_tiles::chip {
namespace {

// Compute and write out a directory for hit rate plots.
void write_chip_hit_rate_directory(TDirectory* parent,
                                   const OccupancyMode& mode,
                                   const ChipSum& chip_sum,
                                   std::uint64_t n_events) {
  const int group = all_groups_index(mode);
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
    const double hit_rate = static_cast<double>(stats.total_hits) / (static_cast<double>(n_events) * kEventWindowSec);
    h_hit_rate->Fill(hit_rate);
  }

  h_hit_rate->Write();
  draw_and_write(dir, h_hit_rate, "c_hit_rate", false, true);
  draw_and_write(dir, chip_sum.h_hits_evt, "c_hits_evt", false, true);
}

// Compute and write out a directory for data rate plots.
void write_chip_data_rate_directory(TDirectory* parent,
                                    const OccupancyMode& mode,
                                    const ChipSum& chip_sum,
                                    std::uint64_t n_events) {
  const int group = all_groups_index(mode);
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
    const double data_rate = stats.total_payload_bits / (static_cast<double>(n_events) * kEventWindowSec * 1.0e9);
    h_data_rate->Fill(data_rate);
  }

  h_data_rate->Write();
  draw_and_write(dir, h_data_rate, "c_data_rate", false, true);
}

}  // namespace

void init_chip_sum(ChipSum& chip_sum, const OccupancyMode& mode) {
  chip_sum.h_hits_evt = new TH1D(
      "h_hits_evt",
      (group_title(mode, all_groups_index(mode)) + ";hits/event/chip;Events").c_str(),
      200,
      0,
      2000);
}

void accumulate_event(ChipSum& chip_sum,
                      const std::array<ChannelHitCountMap, 2>& side_channel_hit_counts,
                      const rates::InsertToLFHCALMapper& mapper) {
  std::unordered_map<rates::VirtualLFHCALPizzaChipID, int, rates::VirtualLFHCALPizzaChipIDHash> chip_active_channel_counts;

  // Loop through both sides.
  for (int side = 0; side < 2; ++side) {

    // Per side, loop thru all channels.
    for (const auto& [channel, count] : side_channel_hit_counts[side]) {
      (void)count;

      // Map the channel onto one of the transverse slices.
      ++chip_active_channel_counts[mapper.pizza_chip(channel, side)];
    }
  }

  // After accumulating all active channels, compute payload bits.
  for (const auto& [chip, active_channel_count] : chip_active_channel_counts) {
    const double payload_bits =
        (kOverheadBits + kBitsPerHit * static_cast<double>(active_channel_count)) * kSamplesPerEvent;

    auto& stats = chip_sum.chips[chip];
    ++stats.total_hits;
    stats.total_payload_bits += payload_bits;
  }

  chip_sum.h_hits_evt->Fill(chip_active_channel_counts.size());
}

void write_output(TDirectory* parent,
                  const OccupancyMode& mode,
                  const ChipSum& chip_sum,
                  std::uint64_t n_events) {
  auto* hit_dir = parent->mkdir("chip_hit_rate");
  auto* data_dir = parent->mkdir("chip_data_rate");
  write_chip_hit_rate_directory(hit_dir, mode, chip_sum, n_events);
  write_chip_data_rate_directory(data_dir, mode, chip_sum, n_events);
}

}  // namespace rates::insert_analysis::lfhcal_tiles::chip
