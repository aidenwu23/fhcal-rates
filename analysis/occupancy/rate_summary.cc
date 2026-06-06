/*

./build/rate_summary -i data/reco -o plots/occupancy/rate_summary.csv

*/

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "classify_hit.h"
#include "utils.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kHitCollection = "LFHCALHits";
constexpr double kDefaultThresholdGeV = 5e-4;
constexpr double kEventWindowSec = 2e-6;

struct Args {
  std::string input_dir;
  std::string output_file;
  double threshold_geV = kDefaultThresholdGeV;
};

struct ColumnStats {
  std::unordered_map<std::uint64_t, std::uint64_t> channel_hits;
  std::uint64_t total_hits = 0;
};

enum class Column : int {
  DIS = 0,
  ElectronBeamBackground = 1,
  ProtonBeamBackground = 2,
  AllSources = 3,
};

constexpr int kNColumns = 4;
const char* kColumnNames[kNColumns] = {
    "DIS",
    "electron_beam_background",
    "proton_beam_background",
    "all_sources",
};

// ----------------------------- handle CLI inputs -----------------------------
void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT_DIR -o OUTPUT.csv [-t THRESHOLD_GEV]\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_dir = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      args.output_file = argv[++i];
    } else if ((arg == "-t" || arg == "--threshold") && i + 1 < argc) {
      args.threshold_geV = std::stod(argv[++i]);
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  if (args.input_dir.empty() || args.output_file.empty()) {
    usage(argv[0]);
    std::exit(1);
  }

  return args;
}

// --------------------------- compute rate numbers ----------------------------
// Highest per-channel rate.
double hottest_channel_rate_hz(const ColumnStats& stats, double total_time_sec) {
  std::uint64_t hottest_hits = 0;
  for (const auto& [channel, hits] : stats.channel_hits) {
    (void)channel;
    hottest_hits = std::max(hottest_hits, hits);
  }
  return total_time_sec > 0.0 ? static_cast<double>(hottest_hits) / total_time_sec : 0.0;
}

// Avg rate per channel.
double mean_channel_rate_hz(const ColumnStats& stats, double total_time_sec) {
  if (stats.channel_hits.empty() || total_time_sec <= 0.0) return 0.0;
  return static_cast<double>(stats.total_hits) /
         (static_cast<double>(stats.channel_hits.size()) * total_time_sec);
}

// Rate summed across all channels.
double total_rate_hz(const ColumnStats& stats, double total_time_sec) {
  return total_time_sec > 0.0 ? static_cast<double>(stats.total_hits) / total_time_sec : 0.0;
}
// -----------------------------------------------------------------------------

// Write output csv.
void write_csv(const fs::path& output_path,
               const std::vector<ColumnStats>& columns,
               std::uint64_t n_events) {
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  std::ofstream out(output_path);
  if (!out) {
    throw std::runtime_error("Failed to open output CSV");
  }

  out << "metric";
  for (int i = 0; i < kNColumns; ++i) out << ',' << kColumnNames[i];
  out << '\n';

  out << std::setprecision(10);

  out << "hottest_channel_avg_hz";
  for (int i = 0; i < kNColumns; ++i) {
    out << ',' << hottest_channel_rate_hz(columns[i], total_time_sec);
  }
  out << '\n';

  out << "channel_avg_hz";
  for (int i = 0; i < kNColumns; ++i) {
    out << ',' << mean_channel_rate_hz(columns[i], total_time_sec);
  }
  out << '\n';

  out << "total_hz";
  for (int i = 0; i < kNColumns; ++i) {
    out << ',' << total_rate_hz(columns[i], total_time_sec);
  }
  out << '\n';
}

}  // namespace

int main(int argc, char* argv[]) {
  const auto args = parse_args(argc, argv);
  const auto files = br::find_root_files(args.input_dir);
  const double threshold_geV = args.threshold_geV;
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_dir << "\n";
    return 1;
  }

  std::vector<ColumnStats> columns(kNColumns);
  std::uint64_t n_events = 0;

  br::FileProgress progress(files.size(), std::cerr);

  // Loop through input files.
  for (const auto& path : files) {
    progress.tick();
    podio::ROOTReader reader;
    reader.openFile(path.string());

    const std::size_t total_events = reader.getEntries("events");

    // Per file, loop through events.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {

      // Grab event + collection.
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!br::has_collection(frame, kHitCollection)) continue;
      ++n_events;

      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop hits.
      for (const auto& hit : hits) {

        // Supress low-energy hits.
        if (hit.getEnergy() <= threshold_geV) continue;

        // Grab dominant contribution and increment stats.
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        const Column source = static_cast<Column>(br::classify_background_class(br::dominant_status(hit)));
        const int source_index = static_cast<int>(source);

        ++columns[source_index].channel_hits[cell_id];
        ++columns[source_index].total_hits;
        ++columns[static_cast<int>(Column::AllSources)].channel_hits[cell_id];
        ++columns[static_cast<int>(Column::AllSources)].total_hits;
      }
    }
  }

  if (n_events == 0) {
    std::cerr << "No events with " << kHitCollection << " found\n";
    return 1;
  }

  try {
    write_csv(args.output_file, columns, n_events);
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << '\n';
    return 1;
  }

  std::cout << "\n";
  return 0;
}
