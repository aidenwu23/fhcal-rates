#pragma once

#include "original_insert/include/family.h"

namespace rates::insert_analysis::original_insert::with_inner_ring {

void init_outputs(Outputs& outputs);
void accumulate_event(Outputs& outputs,
                      const std::vector<rates::insert_analysis::EventHit>& event_hits,
                      const rates::InsertToLFHCALMapper& mapper);
void write_output(TFile& output,
                  const Outputs& outputs,
                  std::uint64_t n_events);

}  // namespace rates::insert_analysis::original_insert::with_inner_ring
