// See shared.cc for actual filtering logic.
#include "lfhcal_tiles/include/no_inner_ring.h"

namespace rates::insert_analysis::lfhcal_tiles::no_inner_ring {

namespace {
OccupancyMode mode() {
  return OccupancyMode{
      "lfhcal_tiles",
      "no_inner_ring",
      true,
  };
}
}  // namespace

void init_outputs(Outputs& outputs) {
  outputs.mode = mode();
  lfhcal_tiles::init_outputs(outputs);
}

void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper) {
  lfhcal_tiles::accumulate_event(outputs, event_hits, mapper);
}

void write_output(TFile& output,
                  const Outputs& outputs,
                  std::uint64_t n_events) {
  lfhcal_tiles::write_output(output, outputs, n_events);
}

}  // namespace rates::insert_analysis::lfhcal_tiles::no_inner_ring
