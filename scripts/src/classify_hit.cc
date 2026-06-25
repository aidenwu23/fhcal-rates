#include "classify_hit.h"

#include <edm4hep/MCParticle.h>

namespace br {

int origin_index(int status) {
  const int offset = status % 1000;
  if (status <= 0 || offset == 0 || (offset >= 5 && offset <= 10)) return -1; // null / reserved
  if (status >= 1    && status < 1000) return 0; // signal / DIS family
  if (status >= 2001 && status < 3000) return 1; // synrad
  if (status >= 3001 && status < 4000) return 2; // eBrem
  if (status >= 4001 && status < 5000) return 3; // eTouschek
  if (status >= 5001 && status < 6000) return 4; // eCoulomb
  if (status >= 6001 && status < 7000) return 5; // pBeamGas
  return 6;                                      // undocumented
}

BackgroundClass classify_background_class(int status) {
  const int offset = status % 1000;
  if (status <= 0 || offset == 0 || (offset >= 5 && offset <= 10)) return BackgroundClass::Other;
  if (status >= 1 && status < 1000) return BackgroundClass::DIS;
  if (status >= 6001 && status < 7000) return BackgroundClass::ProtonBeamBackground;
  if (status >= 2001 && status < 6000) return BackgroundClass::ElectronBeamBackground;
  return BackgroundClass::Other;
}

}  // namespace br

int br::dominant_status(const edm4hep::SimCalorimeterHit& hit) {
  double biggest_contribution = -1.0;
  int best_status = 0;

  // Find dominant status by finding the contribution with the highest energy.
  for (const auto& contribution : hit.getContributions()) {
    const double contribution_energy = contribution.getEnergy();
    if (contribution_energy > biggest_contribution) {
      biggest_contribution = contribution_energy;
      best_status = contribution.getParticle().getGeneratorStatus();
    }
  }

  return best_status;
}
