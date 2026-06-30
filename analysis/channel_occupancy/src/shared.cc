#include "shared.h"

#include <TCanvas.h>
#include <TH2.h>

#include <cmath>
#include <string>

namespace rates::channel_occupancy {
// ----------------------------------------------------------------------------------
// Struct methods
// ----------------------------------------------------------------------------------
double ChannelStats::r() const { return std::hypot(x_mm, y_mm); }

// ----------------------------------------------------------------------------------
// Helpers
// ----------------------------------------------------------------------------------
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

std::vector<double> make_axis_edges(const std::set<double>& coords) {
  if (coords.empty()) return {-kDisplayPaddingMM, kDisplayPaddingMM};

  std::vector<double> values(coords.begin(), coords.end());
  std::vector<double> edges;
  edges.reserve(values.size() + 1);

  if (values.size() == 1) {
    edges.push_back(values.front() - kDisplayPaddingMM);
    edges.push_back(values.front() + kDisplayPaddingMM);
    return edges;
  }

  edges.push_back(values.front() - 0.5 * (values[1] - values[0]));
  for (std::size_t i = 0; i + 1 < values.size(); ++i) {
    edges.push_back(0.5 * (values[i] + values[i + 1]));
  }
  edges.push_back(values.back() + 0.5 * (values.back() - values[values.size() - 2]));
  return edges;
}

AxisEdges2D make_layer_axis_edges(const LayerAccum& layer_accum) {
  std::set<double> x_coords;
  std::set<double> y_coords;
  for (const auto& [channel_id, stats] : layer_accum.channels) {
    (void)channel_id;
    x_coords.insert(stats.x_mm);
    y_coords.insert(stats.y_mm);
  }
  return AxisEdges2D{make_axis_edges(x_coords), make_axis_edges(y_coords)};
}

// ----------------------------------------------------------------------------------
// Initialization
// ----------------------------------------------------------------------------------
void init_layer_accumulations(std::vector<LayerAccum>& layers) {
  layers.assign(kNLayers + 1, {});
  for (int layer = 0; layer <= kNLayers; ++layer) {
    layers[layer].h_hits_evt = new TH1D(
        "h_hits_evt",
        (layer == kAllLayersIndex ? std::string("LFHCAL summed layers;hits/event;Events")
                                  : std::string("LFHCAL layer ") + std::to_string(layer) + ";hits/event;Events").c_str(),
        200, 0, 200);
  }
}

}  // namespace rates::channel_occupancy
