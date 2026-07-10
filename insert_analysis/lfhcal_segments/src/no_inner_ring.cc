// See shared.cc for actual filtering logic.

#include "lfhcal_segments/include/no_inner_ring.h"

namespace rates::insert_analysis::lfhcal_segments::no_inner_ring {

namespace {
OccupancyMode mode() {
  return OccupancyMode{
      "lfhcal_segments",
      "no_inner_ring",
      true,
  };
}
}  // namespace

void init_outputs(Outputs& outputs) {
  outputs.mode = mode();
  lfhcal_segments::init_outputs(outputs);
}

void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper) {
  lfhcal_segments::accumulate_event(outputs, event_hits, mapper);
}

void write_output(TFile& output,
                  const Outputs& outputs,
                  std::uint64_t n_events) {
  lfhcal_segments::write_output(output, outputs, n_events);
}

}  // namespace rates::insert_analysis::lfhcal_segments::no_inner_ring
