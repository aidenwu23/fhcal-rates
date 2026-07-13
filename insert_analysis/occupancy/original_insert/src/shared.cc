#include "original_insert/include/shared.h"

#include <TCanvas.h>
#include <TH2.h>

#include <string>

namespace rates::insert_analysis::original_insert {

int mapped_layer(const OccupancyMode& mode, int physical_layer) {
  // Use zero-based layer indices internally while retaining one-based output labels.
  (void)mode;
  return physical_layer - 1;
}

int group_count(const OccupancyMode& mode) {
  (void)mode;
  return kNLayers;
}

int all_groups_index(const OccupancyMode& mode) {
  return group_count(mode);
}

std::string group_dir_name(const OccupancyMode& mode, int group) {
  if (group == all_groups_index(mode)) {
    return "sum_layers";
  }
  return "layer" + std::to_string(group + 1);
}

std::string group_title(const OccupancyMode& mode, int group) {
  if (group == all_groups_index(mode)) {
    return "Insert all layers";
  }
  return "Insert layer " + std::to_string(group + 1);
}

double channel_mip_energy_gev(const OccupancyMode& mode, int mapped_layer_value) {
  (void)mode;
  (void)mapped_layer_value;
  return kTileMipGeV;
}

bool keep_channel(const OccupancyMode& mode, const OriginalChannelID& channel) {
  if (!mode.drop_inner_ring) return true;

  // Measure the original cell center from the insert's offset inner-ring center.
  const double dx_mm = channel.x_mm - kInnerRingCenterXMM;
  return dx_mm * dx_mm + channel.y_mm * channel.y_mm >= kInnerRingRadiusMM * kInnerRingRadiusMM;
}

std::size_t percentile_index(std::size_t n_values, double percentile) {
  if (n_values == 0) return 0;
  return static_cast<std::size_t>(percentile * static_cast<double>(n_values - 1));
}

// Called when writing one histogram and its canvas to the output file.
void draw_and_write(TDirectory* canvas_dir, TH1* hist, const char* canvas_name, bool logz, bool logy) {
  canvas_dir->cd();
  TCanvas canvas(canvas_name, hist->GetTitle(), 1000, 800);
  if (logz) canvas.SetLogz();
  if (logy) canvas.SetLogy();
  hist->SetStats(false);
  // ROOT uses different draw options for one- and two-dimensional histograms.
  const bool is2d = hist->InheritsFrom(TH2::Class());
  hist->Draw(is2d ? "colz" : "hist");
  canvas.Write();
}

// Called once before the file loop to allocate shared segment products.
void init_segment_sums(std::vector<SegmentSum>& segment_sums, const OccupancyMode& mode) {
  segment_sums.assign(group_count(mode) + 1, {});

  // Each segment gets its own hits-per-event histogram, plus one inclusive slot.
  for (int group = 0; group <= group_count(mode); ++group) {
    segment_sums[group].h_hits_evt = new TH1D(
        "h_hits_evt",
        (group_title(mode, group) + ";hits/event;Events").c_str(),
        200,
        0,
        2000);
  }
}

// Called per event after the hit loop to fill hits-per-event histograms.
void fill_event_histograms(std::vector<SegmentSum>& segment_sums,
                           const OccupancyMode& mode,
                           const std::unordered_map<OriginalChannelID, int, OriginalChannelIDHash>& channel_hit_counts) {
  std::vector<int> hits_per_event(group_count(mode) + 1, 0);

  // Sum this event's hit multiplicity once per segment and once for the inclusive view.
  for (const auto& [channel, count] : channel_hit_counts) {
    const int group = channel.layer;
    if (group < 0 || group >= all_groups_index(mode)) continue;
    hits_per_event[group] += count;
    hits_per_event[all_groups_index(mode)] += count;
  }

  // After the event is summed, fill the corresponding per-segment histograms.
  for (int group = 0; group <= group_count(mode); ++group) {
    segment_sums[group].h_hits_evt->Fill(hits_per_event[group]);
  }
}

}  // namespace rates::insert_analysis::original_insert
