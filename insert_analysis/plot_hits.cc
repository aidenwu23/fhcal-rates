/*

./build/plot_hits -i data/test -o insert_plots/plot_hits.root

*/

#include <TCanvas.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace {
// ----------------------------------------------------------------------------------
// Constants and structs
// ----------------------------------------------------------------------------------

constexpr const char* kCollectionMatch = "HcalEndcapPInsert";

constexpr int kMaxLayersToPlot = 4;
constexpr int kFirstLayer = 1;

// Show a fixed 500 mm radius window around the origin.
constexpr double kWindowRadiusMm = 500.0;
constexpr double kVirtualCellSizeMm = 50.0;
constexpr double kInnerRingCenterXMM = -72.0;
constexpr double kInnerRingRadiusMm = 163.0;
constexpr double kPizzaCenterXMM = 0.5 * (-3.95 - 21.5);
constexpr double kPi = 3.14159265358979323846;
constexpr double kLeftBoundaryAngle = 0.6580536749365931;
constexpr double kRightBoundaryAngle = 0.757476810561768;

const std::array<int, kMaxLayersToPlot> kLayerColors = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};
constexpr int kLeftMarkerStyle = 20;
constexpr int kRightMarkerStyle = 21;
struct Bucket {
  // Unique physical cells and their centers for one layer and insert side.
  std::unordered_set<std::uint64_t> cells;
  std::vector<double> xs;
  std::vector<double> ys;
};

// ----------------------------------------------------------------------------------
// CLI
// ----------------------------------------------------------------------------------
struct Args {
  std::string input_dir;
  std::string output_file;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT_DIR -o OUTPUT.root\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_dir = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      args.output_file = argv[++i];
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input_dir.empty() || args.output_file.empty()) {
    usage(argv[0]);
    std::exit(1);
  }

  return args;
}

// Pick the first insert hit collection that stores SimCalorimeterHit objects.
std::string find_insert_collection(const podio::Frame& frame) {
  for (const auto& name : frame.getAvailableCollections()) {
    if (name.find(kCollectionMatch) != std::string::npos &&
        name.find("Contributions") == std::string::npos) {
      return name;
    }
  }
  return "";
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TH1::AddDirectory(false);

  // Parse arguments and find every ROOT file in the input directory.
  const auto args = parse_args(argc, argv);
  const auto files = rates::find_root_files(args.input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  // Ensure the output directory exists before any event work starts.
  fs::path output_path = args.output_file;
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // The cell ID provides layer and side, while the hit position provides x and y.
  const rates::HcalEndcapPInsertCellIDDecoder decoder;
  std::array<std::array<Bucket, 2>, kMaxLayersToPlot> buckets;

  std::string collection_name;
  rates::FileProgress progress(files.size(), std::cerr);

  // Loop over all files.
  for (const auto& path : files) {
    progress.tick();

    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t total_events = reader.getEntries("events");

    // Loop over all events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));

      // Lock onto the first insert hit collection name and reuse it for later events.
      if (collection_name.empty()) {
        collection_name = find_insert_collection(frame);
      }
      if (collection_name.empty() || !rates::has_collection(frame, collection_name)) continue;

      // Keep one marker per physical cell while retaining its first observed position.
      const auto& event_hits = frame.get<edm4hep::SimCalorimeterHitCollection>(collection_name);

      // Loop over all insert hits.
      for (const auto& hit : event_hits) {
        const std::uint64_t cell_id = static_cast<std::uint64_t>(hit.getCellID());
        const auto cell = decoder.cell(cell_id);
        const int layer_slot = cell.layer - kFirstLayer;
        if (layer_slot < 0 || layer_slot >= kMaxLayersToPlot) continue;
        if (cell.side < 0 || cell.side > 1) continue;

        auto& bucket = buckets[layer_slot][cell.side];
        if (!bucket.cells.insert(cell_id).second) continue;  // Keep one marker per physical cell.

        const auto position = hit.getPosition();
        const double x_mm = static_cast<double>(position.x);
        const double y_mm = static_cast<double>(position.y);

        bucket.xs.push_back(x_mm);
        bucket.ys.push_back(y_mm);
      }
    }
  }

  std::cerr << "\n";

  if (collection_name.empty()) {
    std::cerr << "Failed to find an insert hit collection\n";
    return 1;
  }

  // Print the final unique-cell counts before building the overlay plot.
  std::size_t total_unique_cells = 0;
  for (int layer_slot = 0; layer_slot < kMaxLayersToPlot; ++layer_slot) {
    const int layer = kFirstLayer + layer_slot;
    const std::size_t count = buckets[layer_slot][0].cells.size() + buckets[layer_slot][1].cells.size();
    total_unique_cells += count;
    std::cout << "layer " << layer << ": " << count << " unique cells\n";
  }
  std::cout << "layers " << kFirstLayer << "-" << (kFirstLayer + kMaxLayersToPlot - 1)
            << " total: " << total_unique_cells << " unique cells\n";

  if (total_unique_cells == 0) {
    std::cerr << "No insert hits found in " << args.input_dir << "\n";
    return 1;
  }

  // Open the output ROOT file and create a blank 2D frame for the marker overlay.
  TFile output(args.output_file.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  TH2D frame_hist(
      "h_xy_frame",
      "Insert hits;x [mm];y [mm]",
      100,
      -kWindowRadiusMm,
      kWindowRadiusMm,
      100,
      -kWindowRadiusMm,
      kWindowRadiusMm);
  frame_hist.SetStats(false);

  TCanvas canvas("c_insert_hits_xy", "Insert hits", 1000, 900);
  frame_hist.Draw();

  // Overlay the 5 x 5 cm virtual-channel boundaries used by the rate analysis.
  std::vector<TLine> segmentation_lines;
  segmentation_lines.reserve(2 * static_cast<int>(2.0 * kWindowRadiusMm / kVirtualCellSizeMm + 1.0));
  for (double coordinate = -kWindowRadiusMm; coordinate <= kWindowRadiusMm; coordinate += kVirtualCellSizeMm) {
    segmentation_lines.emplace_back(coordinate, -kWindowRadiusMm, coordinate, kWindowRadiusMm);
    auto& vertical_line = segmentation_lines.back();
    vertical_line.SetLineColor(kGray + 2);
    vertical_line.SetLineStyle(1);
    vertical_line.SetLineWidth(2);
    vertical_line.Draw();

    segmentation_lines.emplace_back(-kWindowRadiusMm, coordinate, kWindowRadiusMm, coordinate);
    auto& horizontal_line = segmentation_lines.back();
    horizontal_line.SetLineColor(kGray + 2);
    horizontal_line.SetLineStyle(1);
    horizontal_line.SetLineWidth(2);
    horizontal_line.Draw();
  }

  const std::array<double, 5> chip_boundary_angles = {
      -kRightBoundaryAngle,
      kRightBoundaryAngle,
      kPi - kLeftBoundaryAngle,
      kPi,
      kPi + kLeftBoundaryAngle,
  };

  std::vector<TLine> chip_lines;
  chip_lines.reserve(chip_boundary_angles.size() + 1);
  chip_lines.emplace_back(kPizzaCenterXMM, -kWindowRadiusMm, kPizzaCenterXMM, kWindowRadiusMm);
  auto& side_line = chip_lines.back();
  side_line.SetLineColor(kRed + 1);
  side_line.SetLineStyle(1);
  side_line.SetLineWidth(2);

  for (const double angle : chip_boundary_angles) {
    const double dx = std::cos(angle);
    const double dy = std::sin(angle);
    double ray_length = (dx > 0.0 ? kWindowRadiusMm - kPizzaCenterXMM
                                  : -kWindowRadiusMm - kPizzaCenterXMM) / dx;
    if (dy != 0.0) {
      const double y_ray_length = (dy > 0.0 ? kWindowRadiusMm : -kWindowRadiusMm) / dy;
      ray_length = std::min(ray_length, y_ray_length);
    }
    chip_lines.emplace_back(
        kPizzaCenterXMM,
        0.0,
        kPizzaCenterXMM + ray_length * dx,
        ray_length * dy);
    auto& line = chip_lines.back();
    line.SetLineColor(kRed + 1);
    line.SetLineStyle(1);
    line.SetLineWidth(2);
  }

  TEllipse inner_ring_cut(
      kInnerRingCenterXMM,
      0.0,
      kInnerRingRadiusMm,
      kInnerRingRadiusMm);
  inner_ring_cut.SetFillStyle(0);
  inner_ring_cut.SetLineColor(kOrange + 7);
  inner_ring_cut.SetLineWidth(2);
  inner_ring_cut.Draw();

  // Keep graphs alive until the canvas is written.
  std::vector<TGraph> graphs;
  graphs.reserve(kMaxLayersToPlot * 2);

  std::array<TGraph*, kMaxLayersToPlot> layer_legend_graphs{};

  // Loop over the requested layers and draw one graph for each side that has cells.
  for (int layer_slot = 0; layer_slot < kMaxLayersToPlot; ++layer_slot) {
    const int layer_value = kFirstLayer + layer_slot;

    // Per layer, loop through left and right side.
    for (int side = 0; side < 2; ++side) {
      const auto& bucket = buckets[layer_slot][side];

      if (bucket.xs.empty()) continue;

      graphs.emplace_back(static_cast<int>(bucket.xs.size()), bucket.xs.data(), bucket.ys.data());
      auto& graph = graphs.back();

      graph.SetName(("g_layer" + std::to_string(layer_value) + "_side" + std::to_string(side)).c_str());
      graph.SetTitle(("Layer " + std::to_string(layer_value)).c_str());
      graph.SetMarkerColor(kLayerColors[layer_slot]);
      graph.SetLineColor(kLayerColors[layer_slot]);
      graph.SetMarkerStyle(side == 0 ? kLeftMarkerStyle : kRightMarkerStyle);
      graph.SetMarkerSize(side == 0 ? 0.8 : 1.0);
      graph.Draw("P SAME");

      if (layer_legend_graphs[layer_slot] == nullptr) layer_legend_graphs[layer_slot] = &graph;
    }
  }

  TLegend legend(0.12, 0.62, 0.32, 0.88);
  legend.SetBorderSize(1);
  legend.SetLineColor(kBlack);
  legend.SetFillColor(kWhite);
  legend.SetFillStyle(1001);

  TGraph left_legend_marker;
  left_legend_marker.SetMarkerColor(kBlack);
  left_legend_marker.SetMarkerStyle(24);
  TGraph right_legend_marker;
  right_legend_marker.SetMarkerColor(kBlack);
  right_legend_marker.SetMarkerStyle(25);
  legend.AddEntry(&left_legend_marker, "left", "p");
  legend.AddEntry(&right_legend_marker, "right", "p");
  legend.AddEntry(static_cast<TObject*>(nullptr), " ", "");
  for (int layer_slot = 0; layer_slot < kMaxLayersToPlot; ++layer_slot) {
    if (layer_legend_graphs[layer_slot] == nullptr) continue;
    const int layer_value = kFirstLayer + layer_slot;
    legend.AddEntry(layer_legend_graphs[layer_slot], ("Layer " + std::to_string(layer_value)).c_str(), "p");
  }

  // Write the frame, marker graphs, and canvases into the output file.
  legend.Draw();

  frame_hist.Write();
  for (auto& graph : graphs) {
    graph.Write();
  }
  canvas.Write();

  TCanvas pizza_canvas("c_insert_hits_xy_pizza", "Insert hits with pizza-chip regions", 1000, 900);
  frame_hist.Draw();
  for (auto& line : segmentation_lines) line.Draw();
  for (auto& line : chip_lines) line.Draw();
  inner_ring_cut.Draw();
  for (auto& graph : graphs) graph.Draw("P SAME");
  legend.Draw();
  pizza_canvas.Write();

  output.Close();
  return 0;
}
