#include "occupancy/include/shared.h"

#include <TCanvas.h>

#include <string>

namespace rates::lfhcal_analysis::occupancy {

std::string layer_name(int layer) {
  // Return the output name for one layer or the inclusive layer group.
  return layer == kAllLayers ? "sum_layers" : "layer" + std::to_string(layer + 1);
}

std::string layer_title(int layer) {
  // Return the plot title for one layer or the inclusive layer group.
  return layer == kAllLayers ? "LFHCAL all readout layers" : "LFHCAL readout layer " + std::to_string(layer + 1);
}

std::size_t percentile_index(std::size_t n_values, double percentile) {
  // Map a percentile fraction onto a sorted vector index.
  if (n_values == 0) return 0;
  return static_cast<std::size_t>(percentile * static_cast<double>(n_values - 1));
}

void draw_and_write(TDirectory* directory, TH1* histogram, const char* canvas_name, bool logx, bool logy) {
  // Draw and write one histogram with the requested axis scales.
  directory->cd();
  TCanvas canvas(canvas_name, histogram->GetTitle(), 1000, 800);
  if (logx) canvas.SetLogx();
  if (logy) canvas.SetLogy();
  histogram->SetStats(false);
  histogram->Draw("hist");
  canvas.Write();
}

}  // namespace rates::lfhcal_analysis::occupancy
