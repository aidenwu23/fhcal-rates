#pragma once

#include <TH1D.h>

#include <string>
#include <vector>

namespace br::origins {

struct LayerEdepHists {
  int layer = -1;
  std::string title;
  std::vector<TH1D*> by_origin;
};

LayerEdepHists make_edep_hists(const std::vector<const char*>& origin_keys,
                               const std::vector<int>& origin_colors,
                               int layer,
                               int bins,
                               const double* edges);

std::vector<TH1D*> make_summed_edep_hists(const std::vector<const char*>& origin_keys,
                                          const std::vector<int>& origin_colors,
                                          int bins,
                                          const double* edges);

void sum_edep_into(std::vector<TH1D*>& summed_hists, const LayerEdepHists& layer_hists);

void fill_edep(LayerEdepHists& hists, int origin_index, double hit_energy);

}  // namespace br::origins
