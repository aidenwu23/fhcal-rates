#include "edep.h"

#include <stdexcept>
#include <string>

namespace rates::origins {
namespace {

// Style hists.
void style(TH1D* hist, int color) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetLineWidth(2);
  hist->SetStats(false);
}

// Make hist title.
std::string edep_title(int layer) {
  if (layer < 0) return "LFHCAL hit origin;E_{dep} [GeV];Hits";
  return "LFHCAL readout layer " + std::to_string(layer) + " hit origin;E_{dep} [GeV];Hits";
}

}  // namespace

// Create the per-origin energy-deposit histograms for one readout layer.
// Called once at startup for each readout layer configuration.
LayerEdepHists make_edep_hists(const std::vector<const char*>& origin_labels,
                               const std::vector<int>& origin_colors,
                               int layer,
                               int bins,
                               const double* edges) {
  if (origin_labels.size() != origin_colors.size()) {
    throw std::invalid_argument("Origin label/color size mismatch");
  }

  LayerEdepHists out;
  out.layer = layer;
  out.title = edep_title(layer);
  out.by_origin.reserve(origin_labels.size());

  const std::string layer_tag = layer < 0 ? "sum" : std::to_string(layer);

  // For each origin (generatorStatus family), build one histogram tagged by layer number.
  for (std::size_t i = 0; i < origin_labels.size(); ++i) {
    const std::string name = "h_edep_hit_layer" + layer_tag + "_" + origin_labels[i];
    auto* hist = new TH1D(name.c_str(), out.title.c_str(), bins, edges);
    style(hist, origin_colors[i]);
    out.by_origin.push_back(hist);
  }

  return out;
}

std::vector<TH1D*> make_summed_edep_hists(const std::vector<const char*>& origin_labels,
                                          const std::vector<int>& origin_colors,
                                          int bins,
                                          const double* edges) {
  // Reuse the same histogram factory for the all-layers-summed view.
  return make_edep_hists(origin_labels, origin_colors, -1, bins, edges).by_origin;
}

// Called after the event loop, once per readout layer.
void sum_edep_into(std::vector<TH1D*>& summed_hists, const LayerEdepHists& layer_hists) {
  // The summed and per-layer containers should hold one histogram per origin.
  if (summed_hists.size() != layer_hists.by_origin.size()) {
    throw std::invalid_argument("Summed/layer histogram size mismatch");
  }

  // Accumulate one layer's contribution into the all-layers histogram.
  for (std::size_t i = 0; i < summed_hists.size(); ++i) {
    summed_hists[i]->Add(layer_hists.by_origin[i]);
  }
}

// Called at the hit level.
void fill_edep(LayerEdepHists& hists, int origin_index, double energy) {
  // Guard against invalid origin indices from the caller.
  if (origin_index < 0 || origin_index >= static_cast<int>(hists.by_origin.size())) {
    throw std::out_of_range("Origin index is out of range");
  }

  if (energy <= 0.0) return;

  hists.by_origin[origin_index]->Fill(energy);
}

}  // namespace rates::origins
