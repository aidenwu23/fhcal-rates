#include "smooth_hists.h"

#include <TH2D.h>

#include <cmath>
#include <string>

namespace br {

void smooth_hist_vertical(TH2D* hist) {
  auto* filled = static_cast<TH2D*>(hist->Clone((std::string(hist->GetName()) + "_filled_vertical").c_str()));
  filled->Reset();

  for (int ix = 2; ix < hist->GetNbinsX(); ++ix) {
    for (int iy = 2; iy < hist->GetNbinsY(); ++iy) {
      if (hist->GetBinContent(ix, iy) > 0.0) continue;

      const double left = hist->GetBinContent(ix - 1, iy);
      const double right = hist->GetBinContent(ix + 1, iy);
      if (left <= 0.0 || right <= 0.0) continue;

      const bool vertical_run =
          hist->GetBinContent(ix, iy - 1) <= 0.0 || hist->GetBinContent(ix, iy + 1) <= 0.0;
      if (!vertical_run) continue;

      double down = 0.0;
      for (int jy = iy - 1; jy >= 1; --jy) {
        down = hist->GetBinContent(ix, jy);
        if (down > 0.0) break;
      }

      double up = 0.0;
      for (int jy = iy + 1; jy <= hist->GetNbinsY(); ++jy) {
        up = hist->GetBinContent(ix, jy);
        if (up > 0.0) break;
      }

      if (down <= 0.0 || up <= 0.0) continue;

      const double log_mean = std::exp(0.25 * (std::log(left) + std::log(right) + std::log(down) + std::log(up)));
      filled->SetBinContent(ix, iy, log_mean);
    }
  }

  hist->Add(filled, 1.0);
  delete filled;
}

void smooth_hist_neighbhors(TH2D* hist) {
  auto* filled = static_cast<TH2D*>(hist->Clone((std::string(hist->GetName()) + "_filled_neighbors").c_str()));
  filled->Reset();

  for (int ix = 2; ix < hist->GetNbinsX(); ++ix) {
    for (int iy = 2; iy < hist->GetNbinsY(); ++iy) {
      if (hist->GetBinContent(ix, iy) > 0.0) continue;

      const double left = hist->GetBinContent(ix - 1, iy);
      const double right = hist->GetBinContent(ix + 1, iy);
      const double down = hist->GetBinContent(ix, iy - 1);
      const double up = hist->GetBinContent(ix, iy + 1);
      if (left <= 0.0 || right <= 0.0 || down <= 0.0 || up <= 0.0) continue;

      const double log_mean = std::exp(0.25 * (std::log(left) + std::log(right) + std::log(down) + std::log(up)));
      filled->SetBinContent(ix, iy, log_mean);
    }
  }

  hist->Add(filled, 1.0);
  delete filled;
}

void smooth_hist_horizontal(TH2D* hist) {
  auto* filled = static_cast<TH2D*>(hist->Clone((std::string(hist->GetName()) + "_filled_horizontal").c_str()));
  filled->Reset();

  for (int ix = 2; ix < hist->GetNbinsX(); ++ix) {
    for (int iy = 2; iy < hist->GetNbinsY(); ++iy) {
      if (hist->GetBinContent(ix, iy) > 0.0) continue;

      const double down = hist->GetBinContent(ix, iy - 1);
      const double up = hist->GetBinContent(ix, iy + 1);
      if (down <= 0.0 || up <= 0.0) continue;

      const bool horizontal_run =
          hist->GetBinContent(ix - 1, iy) <= 0.0 || hist->GetBinContent(ix + 1, iy) <= 0.0;
      if (!horizontal_run) continue;

      double left = 0.0;
      for (int jx = ix - 1; jx >= 1; --jx) {
        left = hist->GetBinContent(jx, iy);
        if (left > 0.0) break;
      }

      double right = 0.0;
      for (int jx = ix + 1; jx <= hist->GetNbinsX(); ++jx) {
        right = hist->GetBinContent(jx, iy);
        if (right > 0.0) break;
      }

      if (left <= 0.0 || right <= 0.0) continue;

      const double log_mean = std::exp(0.25 * (std::log(left) + std::log(right) + std::log(down) + std::log(up)));
      filled->SetBinContent(ix, iy, log_mean);
    }
  }

  hist->Add(filled, 1.0);
  delete filled;
}

}  // namespace br
