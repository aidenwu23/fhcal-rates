#pragma once

#include <TH1D.h>

#include <edm4hep/Vector3f.h>

#include <string>
#include <vector>

namespace rates::origins {

struct ThresholdEtaHists {
  double threshold_geV = 0.0;
  std::string directory;
  std::string canvas_name;
  std::string title;
  std::vector<TH1D*> by_origin;
};

double eta(const edm4hep::Vector3f& position);

ThresholdEtaHists make_eta_hists(const std::vector<const char*>& origin_labels,
                                 const std::vector<int>& origin_colors,
                                 double threshold_geV,
                                 int bins = 240,
                                 double eta_min = -12.0,
                                 double eta_max = 12.0);

void fill_eta(ThresholdEtaHists& hists,
              int origin_index,
              const edm4hep::Vector3f& position,
              double hit_energy);

}  // namespace rates::origins
