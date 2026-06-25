#include "classify_hit.h"

#include <edm4hep/MCParticle.h>

namespace br {

int origin_index(int status) {
  // Even when generatorStatus is a null-like entry such as 0, 2000, 3000, ...
  // we still keep it in its broad source family because those shower byproducts
  // can contribute to detector rates and occupancy.
  if (status >= 0    && status < 1000) return 0; // signal / DIS family
  if (status >= 2000 && status < 3000) return 1; // synrad
  if (status >= 3000 && status < 4000) return 2; // eBrem
  if (status >= 4000 && status < 5000) return 3; // eTouschek
  if (status >= 5000 && status < 6000) return 4; // eCoulomb
  if (status >= 6000 && status < 7000) return 5; // pBeamGas
  return 6;                                      // undocumented

  /*
  // More specific status-kind version:
  // treat offset 0 and offsets 5-10 within each shifted block as null / reserved,
  // and exclude them from the broad-family origin assignment.
  const int offset = status % 1000;
  if (status <= 0 || offset == 0 || (offset >= 5 && offset <= 10)) return -1;
  if (status >= 1    && status < 1000) return 0;
  if (status >= 2001 && status < 3000) return 1;
  if (status >= 3001 && status < 4000) return 2;
  if (status >= 4001 && status < 5000) return 3;
  if (status >= 5001 && status < 6000) return 4;
  if (status >= 6001 && status < 7000) return 5;
  return 6;
  */
}

BackgroundClass classify_background_class(int status) {
  if (status >= 0 && status < 1000) return BackgroundClass::DIS;
  if (status >= 6000 && status < 7000) return BackgroundClass::ProtonBeamBackground;
  if (status >= 2000 && status < 6000) return BackgroundClass::ElectronBeamBackground;
  return BackgroundClass::Other;

  /*
  // More specific status-kind version:
  // offset 0 and offsets 5-10 are treated as null / reserved and mapped to Other.
  const int offset = status % 1000;
  if (status <= 0 || offset == 0 || (offset >= 5 && offset <= 10)) return BackgroundClass::Other;
  if (status >= 1 && status < 1000) return BackgroundClass::DIS;
  if (status >= 6001 && status < 7000) return BackgroundClass::ProtonBeamBackground;
  if (status >= 2001 && status < 6000) return BackgroundClass::ElectronBeamBackground;
  return BackgroundClass::Other;
  */
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
