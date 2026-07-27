#pragma once

#include <TGraph.h>
#include <TH1.h>
#include <TMultiGraph.h>

#include <podio/Frame.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <vector>

namespace rates {

bool has_collection(const podio::Frame& frame, const std::string& name);
std::vector<std::filesystem::path> find_root_files(const std::string& input_dir);
std::vector<double> log_edges(int bins, double low, double high);
double mip_energy_gev(int layer);

// Pad plotted axes beyond the largest observed coordinates.
inline void pad_axes(TGraph& graph) {
  if (graph.GetN() == 0) return;

  double max_x = graph.GetPointX(0);
  double max_y = graph.GetPointY(0);
  for (int point = 1; point < graph.GetN(); ++point) {
    if (graph.GetPointX(point) > max_x) max_x = graph.GetPointX(point);
    if (graph.GetPointY(point) > max_y) max_y = graph.GetPointY(point);
  }

  graph.GetXaxis()->SetLimits(0.0, 1.1 * max_x);
  graph.SetMaximum(1.1 * max_y);
}

inline void pad_axes(TMultiGraph& graphs) {
  double max_x = 0.0;
  double max_y = 0.0;
  bool has_points = false;

  for (int index = 0; index < graphs.GetListOfGraphs()->GetSize(); ++index) {
    auto* graph = static_cast<TGraph*>(graphs.GetListOfGraphs()->At(index));
    for (int point = 0; point < graph->GetN(); ++point) {
      if (!has_points || graph->GetPointX(point) > max_x) max_x = graph->GetPointX(point);
      if (!has_points || graph->GetPointY(point) > max_y) max_y = graph->GetPointY(point);
      has_points = true;
    }
  }
  if (!has_points) return;

  graphs.GetXaxis()->SetLimits(0.0, 1.1 * max_x);
  graphs.SetMaximum(1.1 * max_y);
}

inline void pad_axes(TH1& histogram) {
  int max_x_bin = 0;
  int max_y_bin = 0;

  for (int x_bin = 1; x_bin <= histogram.GetNbinsX(); ++x_bin) {
    for (int y_bin = 1; y_bin <= std::max(histogram.GetNbinsY(), 1); ++y_bin) {
      if (histogram.GetBinContent(histogram.GetBin(x_bin, y_bin)) == 0.0) continue;
      max_x_bin = x_bin;
      max_y_bin = std::max(max_y_bin, y_bin);
    }
  }
  if (max_x_bin == 0) return;

  histogram.GetXaxis()->SetRangeUser(histogram.GetXaxis()->GetXmin(),
                                     1.1 * histogram.GetXaxis()->GetBinUpEdge(max_x_bin));
  if (histogram.GetDimension() == 1) {
    histogram.SetMaximum(1.1 * histogram.GetMaximum());
  } else {
    histogram.GetYaxis()->SetRangeUser(histogram.GetYaxis()->GetXmin(),
                                       1.1 * histogram.GetYaxis()->GetBinUpEdge(max_y_bin));
  }
}

template <typename Histograms>
void pad_axes(TH1& axes, const Histograms& histograms) {
  int max_x_bin = 0;
  int max_y_bin = 0;
  double max_content = 0.0;

  for (const auto* histogram : histograms) {
    for (int x_bin = 1; x_bin <= histogram->GetNbinsX(); ++x_bin) {
      for (int y_bin = 1; y_bin <= std::max(histogram->GetNbinsY(), 1); ++y_bin) {
        const double content = histogram->GetBinContent(histogram->GetBin(x_bin, y_bin));
        if (content == 0.0) continue;
        max_x_bin = std::max(max_x_bin, x_bin);
        max_y_bin = std::max(max_y_bin, y_bin);
        max_content = std::max(max_content, content);
      }
    }
  }
  if (max_x_bin == 0) return;

  axes.GetXaxis()->SetRangeUser(axes.GetXaxis()->GetXmin(),
                                1.1 * axes.GetXaxis()->GetBinUpEdge(max_x_bin));
  if (axes.GetDimension() == 1) {
    axes.SetMaximum(1.1 * max_content);
  } else {
    axes.GetYaxis()->SetRangeUser(axes.GetYaxis()->GetXmin(),
                                  1.1 * axes.GetYaxis()->GetBinUpEdge(max_y_bin));
  }
}

class FileProgress {
 public:
  FileProgress(std::size_t total_files, std::ostream& os);
  void tick();

 private:
  std::size_t total_files_ = 0;
  std::size_t current_file_ = 0;
  std::ostream& os_;
};

}  // namespace rates
