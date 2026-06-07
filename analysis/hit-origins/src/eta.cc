#include "eta.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace br::origins {
namespace {

// Style hists.
void style(TH1D* hist, int color) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetLineWidth(2);
  hist->SetStats(false);
}

// Convert a threshold value into a tag for names/directories.
std::string threshold_tag(double threshold_geV) {
  if (threshold_geV == 0.0) return "all";

  std::ostringstream out;
  out << std::scientific << std::setprecision(0) << threshold_geV;
  auto tag = out.str();
  for (char& c : tag) {
    if (c == '+') c = 'p';
    if (c == '-') c = 'm';
    if (c == '.') c = 'p';
  }
  return tag;
}

// Build the title that describes which energy threshold was applied.
std::string eta_title(double threshold_geV) {
  if (threshold_geV == 0.0) return "LFHCAL hit origin;#eta;Hits";

  std::ostringstream out;
  out << "LFHCAL hit origin (E_{dep} > " << std::scientific << std::setprecision(1)
      << threshold_geV << " GeV);#eta;Hits";
  return out.str();
}

}  // namespace

// Compute pseudorapidity from the hit position.
double eta(const edm4hep::Vector3f& position) {
  const double pt = std::hypot(position.x, position.y);
  if (pt <= 0.0) return 0.0;
  return std::asinh(position.z / pt);
}

// Called once at startup for each threshold configuration.
ThresholdEtaHists make_eta_hists(const std::vector<const char*>& origin_labels,
                                 const std::vector<int>& origin_colors,
                                 double threshold_geV,
                                 int bins,
                                 double eta_min,
                                 double eta_max) {
  if (origin_labels.size() != origin_colors.size()) {
    throw std::invalid_argument("Origin label/color size mismatch");
  }

  ThresholdEtaHists out;
  out.threshold_geV = threshold_geV;
  const auto tag = threshold_tag(threshold_geV);
  out.directory = threshold_geV == 0.0 ? "Eta_hit" : "Eta_hit_" + tag;
  out.canvas_name = threshold_geV == 0.0 ? "c_eta_hit" : "c_eta_hit_" + tag;
  out.title = eta_title(threshold_geV);
  out.by_origin.reserve(origin_labels.size());

  // For each origin, build one eta histogram for this threshold choice.
  for (std::size_t i = 0; i < origin_labels.size(); ++i) {
    const std::string name = "h_eta_hit_" + tag + "_" + origin_labels[i];
    auto* hist = new TH1D(name.c_str(), out.title.c_str(), bins, eta_min, eta_max);
    style(hist, origin_colors[i]);
    out.by_origin.push_back(hist);
  }

  return out;
}

// Called once per hit for each eta-threshold configuration.
void fill_eta(ThresholdEtaHists& hists,
              int origin_index,
              const edm4hep::Vector3f& position,
              double hit_energy) {
  // Guard against invalid origin indices from the caller.
  if (origin_index < 0 || origin_index >= static_cast<int>(hists.by_origin.size())) {
    throw std::out_of_range("Origin index is out of range");
  }

  // Skip hits that fail this threshold configuration.
  if (hit_energy <= hists.threshold_geV) return;
  hists.by_origin[origin_index]->Fill(eta(position));
}

}  // namespace br::origins
