#include "classify_hit.h"

#include <edm4hep/MCParticle.h>

namespace br {

int origin_index(int status) {
  if (status >= 2000 && status < 3000) return 1;
  if (status >= 3000 && status < 4000) return 2;
  if (status >= 4000 && status < 5000) return 3;
  if (status >= 5000 && status < 6000) return 4;
  if (status >= 6000 && status < 7000) return 5;
  return 0;
}

BackgroundClass classify_background_class(int status) {
  if (status >= 6000 && status < 7000) return BackgroundClass::ProtonBeamBackground;
  if (status >= 2000 && status < 6000) return BackgroundClass::ElectronBeamBackground;
  return BackgroundClass::DIS;
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
