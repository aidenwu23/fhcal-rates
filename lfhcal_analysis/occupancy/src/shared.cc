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

int channel_type(int layer) {
  return layer < 2 ? 0 : 1;
}

std::string channel_type_name(int type) {
  return type == 0 ? "5_tile" : "10_tile";
}

std::string channel_type_title(int type) {
  return type == 0 ? "LFHCAL 5-tile channels" : "LFHCAL 10-tile channels";
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
  histogram->SetLineWidth(2);
  histogram->Draw("hist");
  rates::pad_axes(*histogram);
  canvas.Write();
}

}  // namespace rates::lfhcal_analysis::occupancy
