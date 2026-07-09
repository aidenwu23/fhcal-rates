#include "shared.h"

#include <TCanvas.h>
#include <TH2.h>

#include <array>
#include <string>

namespace rates::insert_occupancy {

// Called when mapping a physical layer number into one longitudinal segment.
int segment_index(int layer) {
  // Each segment covers a fixed inclusive layer range.
  for (int segment = 0; segment < kNSegments; ++segment) {
    if (layer >= kSegmentFirstLayers[segment] && layer <= kSegmentLastLayers[segment]) return segment;
  }
  return -1;
}

// Called when naming one longitudinal segment directory.
const char* segment_dir_name(int segment) {
  static constexpr std::array<const char*, kNSegments + 1> kNames = {
      "segment1",
      "segment2",
      "segment3",
      "segment4",
      "segment5",
      "segment6",
      "segment7",
      "sum_segments",
  };
  return kNames[segment];
}

// Called when building one longitudinal segment title.
std::string segment_title(int segment) {
  if (segment == kAllSegmentsIndex) return "Summed segments";
  return "Layers " + std::to_string(kSegmentFirstLayers[segment]) + "-" + std::to_string(kSegmentLastLayers[segment]);
}

// Called when getting the number of physical layers inside one longitudinal segment.
int segment_nlayers(int segment) {
  if (segment < 0 || segment >= kNSegments) return 0;
  return kSegmentLastLayers[segment] - kSegmentFirstLayers[segment] + 1;
}

double channel_mip_energy_gev(int segment) {
  return static_cast<double>(segment_nlayers(segment)) * kTileMipGeV;
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
void init_segment_sums(std::vector<SegmentSum>& segment_sums) {
  segment_sums.assign(kNSegments + 1, {});

  // Each segment gets its own hits-per-event histogram, plus one inclusive slot.
  for (int segment = 0; segment <= kNSegments; ++segment) {
    segment_sums[segment].h_hits_evt = new TH1D(
        "h_hits_evt",
        (segment_title(segment) + ";hits/event;Events").c_str(),
        200,
        0,
        2000);
  }
}

// Called per event after the hit loop to fill hits-per-event histograms.
void fill_event_histograms(std::vector<SegmentSum>& segment_sums,
                           const std::unordered_map<rates::VirtualLFHCALChannelID, int, rates::VirtualLFHCALChannelIDHash>& channel_hit_counts) {
  std::array<int, kNSegments + 1> hits_per_event{};

  // Sum this event's hit multiplicity once per segment and once for the inclusive view.
  for (const auto& [channel, count] : channel_hit_counts) {
    const int segment = segment_index(channel.layer);
    if (segment < 0) continue;
    hits_per_event[segment] += count;
    hits_per_event[kAllSegmentsIndex] += count;
  }

  // After the event is summed, fill the corresponding per-segment histograms.
  for (int segment = 0; segment <= kNSegments; ++segment) {
    segment_sums[segment].h_hits_evt->Fill(hits_per_event[segment]);
  }
}

}  // namespace rates::insert_occupancy
