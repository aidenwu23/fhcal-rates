// See shared.cc for actual filtering logic.
#include "physical_layers/include/with_inner_ring.h"

namespace rates::insert_analysis::physical_layers::with_inner_ring {

namespace {
OccupancyMode mode() {
  return OccupancyMode{
      "physical_layers",
      "with_inner_ring",
      false,
  };
}
}  // namespace

void init_outputs(Outputs& outputs) {
  outputs.mode = mode();
  physical_layers::init_outputs(outputs);
}

void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper) {
  physical_layers::accumulate_event(outputs, event_hits, mapper);
}

void write_output(TFile& output,
                  const Outputs& outputs,
                  std::uint64_t n_events) {
  physical_layers::write_output(output, outputs, n_events);
}

}  // namespace rates::insert_analysis::physical_layers::with_inner_ring
