// See shared.cc for actual filtering logic.

#include "full_lfhcal/include/no_inner_ring.h"

namespace rates::insert_analysis::full_lfhcal::no_inner_ring {

namespace {
OccupancyMode mode() {
  return OccupancyMode{
      "full_lfhcal",
      "no_inner_ring",
      true,
  };
}
}  // namespace

void init_outputs(Outputs& outputs) {
  outputs.mode = mode();
  full_lfhcal::init_outputs(outputs);
}

void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper) {
  full_lfhcal::accumulate_event(outputs, event_hits, mapper);
}

void write_output(TFile& output,
                  const Outputs& outputs,
                  std::uint64_t n_events) {
  full_lfhcal::write_output(output, outputs, n_events);
}

}  // namespace rates::insert_analysis::full_lfhcal::no_inner_ring
