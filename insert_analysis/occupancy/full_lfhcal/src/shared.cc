#include "full_lfhcal/include/shared.h"

#include <TCanvas.h>
#include <TH2.h>

#include <array>
#include <string>

namespace rates::insert_analysis::full_lfhcal {

// Called when mapping a physical layer number into one longitudinal segment.
int segment_index(int layer) {
  // Each segment covers a fixed inclusive layer range.
  for (int segment = 0; segment < kNSegments; ++segment) {
    if (layer >= kSegmentFirstLayers[segment] && layer <= kSegmentLastLayers[segment]) return segment;
  }
  return -1;
}

int mapped_layer(const OccupancyMode& mode, int physical_layer) {
  // Store the zero-based segment index as the shared virtual-channel layer.
  (void)mode;
  const int segment = segment_index(physical_layer);
  return segment;
}

int group_count(const OccupancyMode& mode) {
  (void)mode;
  return kNSegments;
}

int all_groups_index(const OccupancyMode& mode) {
  return group_count(mode);
}

std::string group_dir_name(const OccupancyMode& mode, int group) {
  if (group == all_groups_index(mode)) {
    return "sum_segments";
  }
  return "segment" + std::to_string(group + 1);
}

std::string group_title(const OccupancyMode& mode, int group) {
  if (group == all_groups_index(mode)) {
    return "Insert all segments (full LFHCal readout)";
  }
  return "Insert layers " + std::to_string(kSegmentFirstLayers[group]) + "-" + std::to_string(kSegmentLastLayers[group]) +
         " (full LFHCal readout)";
}

double channel_mip_energy_gev(const OccupancyMode& mode, int mapped_layer_value) {
  (void)mode;
  if (mapped_layer_value < 0 || mapped_layer_value >= kNSegments) return 0.0;
  return static_cast<double>(kSegmentLastLayers[mapped_layer_value] -
                             kSegmentFirstLayers[mapped_layer_value] + 1) *
         kTileMipGeV;
}

bool keep_channel(const OccupancyMode& mode, const rates::VirtualLFHCALChannelID& channel) {
  if (!mode.drop_inner_ring) return true;

  // Measure the virtual-channel center from the insert's offset inner-ring center.
  const double x_mm = (static_cast<double>(channel.ix) + 0.5) * kVirtualCellSizeMM;
  const double y_mm = (static_cast<double>(channel.iy) + 0.5) * kVirtualCellSizeMM;
  const double dx_mm = x_mm - kInnerRingCenterXMM;
  return dx_mm * dx_mm + y_mm * y_mm >= kInnerRingRadiusMM * kInnerRingRadiusMM;
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

// Called when turning filled x or y coordinates into bin edges.
std::vector<double> make_axis_edges(const std::set<double>& coords) {
  // Fall back to one nominal bin when this map never received any channels.
  if (coords.empty()) return {-0.5 * kVirtualCellSizeMM, 0.5 * kVirtualCellSizeMM};

  std::vector<double> values(coords.begin(), coords.end());
  std::vector<double> edges;
  edges.reserve(values.size() + 1);

  // Use the stored channel centers to build one cell-width bin around each coordinate.
  edges.push_back(values.front() - 0.5 * kVirtualCellSizeMM);
  for (double value : values) {
    edges.push_back(value + 0.5 * kVirtualCellSizeMM);
  }
  return edges;
}

// Called when building x and y bin edges for one segment map.
AxisEdges2D make_segment_axis_edges(const SegmentSum& segment_sum) {
  std::set<double> x_coords;
  std::set<double> y_coords;

  // Collect unique channel centers so the output map follows the occupied geometry.
  for (const auto& [channel_id, stats] : segment_sum.channels) {
    (void)channel_id;
    x_coords.insert(stats.x_mm);
    y_coords.insert(stats.y_mm);
  }
  return AxisEdges2D{make_axis_edges(x_coords), make_axis_edges(y_coords)};
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
                           const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& channel_hit_counts) {
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

}  // namespace rates::insert_analysis::full_lfhcal
