#pragma once

#include <edm4hep/SimCalorimeterHit.h>

namespace br {

enum class BackgroundClass : int {
  DIS = 0,
  ElectronBeamBackground = 1,
  ProtonBeamBackground = 2,
  Other = 3,
};

int origin_index(int status);
BackgroundClass classify_background_class(int status);
int dominant_status(const edm4hep::SimCalorimeterHit& hit);

}  // namespace br
