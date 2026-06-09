#pragma once

class TH2D;

namespace br {

void smooth_hist_vertical(TH2D* hist);
void smooth_hist_neighbhors(TH2D* hist);
void smooth_hist_horizontal(TH2D* hist);

}  // namespace br
