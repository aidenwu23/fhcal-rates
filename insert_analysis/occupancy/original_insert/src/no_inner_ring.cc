// See shared.cc for actual filtering logic.
#include "original_insert/include/no_inner_ring.h"

namespace rates::insert_analysis::original_insert::no_inner_ring {

namespace {
OccupancyMode mode() {
  return OccupancyMode{
      "original_insert",
      "no_inner_ring",
      true,
  };
}
}  // namespace

void init_outputs(Outputs& outputs) {
  outputs.mode = mode();
  original_insert::init_outputs(outputs);
}

void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits) {
  original_insert::accumulate_event(outputs, event_hits);
}

void write_output(TFile& output,
                  const Outputs& outputs,
                  std::uint64_t n_events) {
  original_insert::write_output(output, outputs, n_events);
}

}  // namespace rates::insert_analysis::original_insert::no_inner_ring
