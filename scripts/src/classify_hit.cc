#include "classify_hit.h"

#include <edm4hep/MCParticle.h>

namespace br {

int origin_index(int status) {
  if (status >= 0    && status < 1000) return 0; // signal / DIS family
  if (status >= 2000 && status < 3000) return 1; // synrad
  if (status >= 3000 && status < 4000) return 2; // eBrem
  if (status >= 4000 && status < 5000) return 3; // eTouschek
  if (status >= 5000 && status < 6000) return 4; // eCoulomb
  if (status >= 6000 && status < 7000) return 5; // pBeamGas
  return 6;                                      // other
}

BackgroundClass classify_background_class(int status) {
  if (status == 1 || status == 2) return BackgroundClass::DIS;
  if (status >= 6000 && status < 7000) return BackgroundClass::ProtonBeamBackground;
  if (status >= 2000 && status < 6000) return BackgroundClass::ElectronBeamBackground;
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
