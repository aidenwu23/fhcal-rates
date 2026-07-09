#include "chip.h"

#include "shared.h"

#include <TH1D.h>

#include <array>
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
// Called when writing one layer of virtual-chip rate histograms.
void write_layer_directory(TDirectory* parent,
                           const LayerAccum& layer_accum,
                           int layer,
                           std::uint64_t n_events) {
  const std::string dir_name = segment_dir_name(layer);
  auto* dir = parent->mkdir(dir_name.c_str());
  dir->cd();

  auto* hist_dir = dir->mkdir("hists");
  hist_dir->cd();
  layer_accum.h_hits_evt->Write();

  auto* h_hit_rate = new TH1D(
      "h_hit_rate",
      (segment_title(layer) + ";rate [Hz/virtual chip];virtual chips").c_str(),
      200,
      0.0,
      2.0e5);

  auto* h_data_rate = new TH1D(
      "h_data_rate",
      (segment_title(layer) + ";tail data rate [Gb/s/virtual chip];virtual chips").c_str(),
      200,
      0.0,
      0.05);

  // Convert each accumulated virtual chip into one full-sample rate entry.
  for (const auto& [chip, stats] : layer_accum.chips) {
    (void)chip;
    const double hit_rate = static_cast<double>(stats.total_hits) /
                            (static_cast<double>(n_events) * rates::insert_occupancy::kEventWindowSec);
    const double data_rate = stats.total_bits /
                             (static_cast<double>(n_events) * rates::insert_occupancy::kEventWindowSec * 1.0e9);
    h_hit_rate->Fill(hit_rate);
    h_data_rate->Fill(data_rate);
  }

  h_hit_rate->Write();
  h_data_rate->Write();

  rates::insert_occupancy::draw_and_write(dir, h_hit_rate, "c_hit_rate", false, true);
  rates::insert_occupancy::draw_and_write(dir, h_data_rate, "c_data_rate", false, true);
  rates::insert_occupancy::draw_and_write(dir, layer_accum.h_hits_evt, "c_hits_evt", false, true);
}

}  // namespace

// Called once before the file loop to allocate chip-level layer products.
void init_layer_accumulations(std::vector<LayerAccum>& layers) {
  layers.assign(rates::insert_occupancy::kNSegments + 1, {});

  // Each segment gets its own chips-per-event histogram, plus one inclusive slot.
  for (int layer = 0; layer <= rates::insert_occupancy::kNSegments; ++layer) {
    layers[layer].h_hits_evt = new TH1D(
        "h_hits_evt",
        (segment_title(layer) + ";hits/event/chip;Events").c_str(),
        200,
        0,
        2000);
  }
}

// Called per event after the hit loop to accumulate chip-level occupancy.
void accumulate_event(std::vector<LayerAccum>& layers,
                      const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper) {
  std::unordered_map<rates::VirtualLFHCALChipID, int, rates::VirtualLFHCALChipIDHash> chip_active_channels;
  std::array<int, rates::insert_occupancy::kNSegments + 1> hits_per_event{};

  // First collapse the event's fired virtual channels into per-chip active-channel counts.
  for (const auto& [channel, count] : event_hits) {
    (void)count;
    const auto chip = mapper.chip(channel);
    ++chip_active_channels[chip];
  }

  // Then convert each fired chip into occupancy and data-volume counters.
  for (const auto& [chip, active_channel_count] : chip_active_channels) {
    const double event_bits =
        (rates::insert_occupancy::kOverheadBits +
         rates::insert_occupancy::kBitsPerHit * static_cast<double>(active_channel_count)) *
        rates::insert_occupancy::kSamplesPerEvent;

    const int segment = rates::insert_occupancy::kAllSegmentsIndex;

    auto& layer_accum = layers[segment];
    auto& stats = layer_accum.chips[chip];
    stats.x_mm = (static_cast<double>(chip.ix) + 0.5) * kChipSizeMM;
    stats.y_mm = (static_cast<double>(chip.iy) + 0.5) * kChipSizeMM;
    ++stats.total_hits;         // Count one fired event for this virtual chip.
    stats.total_bits += event_bits;  // Accumulate the chip payload from this event.

    ++hits_per_event[segment];
  }

  // Fill the per-event chip multiplicity histograms once the event is fully summed.
  for (int layer = 0; layer <= rates::insert_occupancy::kNSegments; ++layer) {
    layers[layer].h_hits_evt->Fill(hits_per_event[layer]);
  }
}

// Called once after the event loop to write virtual-chip products.
void write_output(TFile& output,
                  const std::vector<LayerAccum>& layers,
                  std::uint64_t n_events) {
  auto* parent = output.mkdir("chip");

  // Write one directory per segment plus one inclusive summed-segment directory.
  for (int layer = 0; layer <= rates::insert_occupancy::kAllSegmentsIndex; ++layer) {
    write_layer_directory(parent, layers[layer], layer, n_events);
  }
}

}  // namespace rates::insert_occupancy::chip
