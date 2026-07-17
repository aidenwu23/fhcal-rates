/*

root -l -b -q 'insert_analysis/overhead_estimate.C("data/test")'

*/

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kThresholdGeV = 0.5 * 4e-4;
constexpr double kEventWindowSec = 2e-6;
constexpr double kVirtualCellSizeMM = 50.0;
constexpr double kHoleCenterXMM = -72.0;
constexpr double kOriginalVetoRadiusMM = 163.0;
constexpr double kLFHCalVetoRadiusMM = 163.0;
constexpr std::uint64_t kChannelsPerChip = 36;
constexpr double kOverheadBits = 128.0;
constexpr double kBitsPerHit = 32.0;
constexpr double kSamplesPerEvent = 4.0;

struct VirtualChannel {
  int layer = 0;
  int ix = 0;
  int iy = 0;

  bool operator==(const VirtualChannel& other) const {
    return layer == other.layer && ix == other.ix && iy == other.iy;
  }
};

struct VirtualChannelHash {
  std::size_t operator()(const VirtualChannel& channel) const {
    std::size_t h = 0;
    const auto mix = [&](int value) {
      h ^= std::hash<int>{}(value) + 0x9e3779b9 + (h << 6) + (h >> 2);
    };
    mix(channel.layer);
    mix(channel.ix);
    mix(channel.iy);
    return h;
  }
};

struct Estimate {
  std::uint64_t active_channels = 0;
  std::uint64_t packed_chips = 0;
  double payload_bits = 0.0;
  double packed_bits = 0.0;
  double per_channel_overhead_bits = 0.0;
};

std::vector<std::string> find_root_files(const std::string& input_dir) {
  std::vector<std::string> files;
  for (const auto& entry : fs::recursive_directory_iterator(input_dir)) {
    if (!entry.is_regular_file()) continue;
    if (entry.path().extension() == ".root") files.push_back(entry.path().string());
  }
  std::sort(files.begin(), files.end());
  return files;
}

std::string find_insert_collection(const podio::Frame& frame) {
  for (const auto& name : frame.getAvailableCollections()) {
    if (name.find("HcalEndcapPInsert") != std::string::npos &&
        name.find("Contributions") == std::string::npos) {
      return name;
    }
  }
  return "";
}

int decode_side(std::uint64_t cell_id) {
  return static_cast<int>((cell_id >> 8) & 0x1);
}

int decode_layer(std::uint64_t cell_id) {
  return static_cast<int>((cell_id >> 9) & 0xff);
}

bool inside_veto(double x_mm, double y_mm, double radius_mm) {
  const double dx_mm = x_mm - kHoleCenterXMM;
  return dx_mm * dx_mm + y_mm * y_mm < radius_mm * radius_mm;
}

void add_event(Estimate& estimate, const std::array<std::uint64_t, 2>& active_by_side) {
  std::uint64_t active_channels = 0;
  std::uint64_t packed_chips = 0;
  for (int side = 0; side < 2; ++side) {
    active_channels += active_by_side[side];
    packed_chips += (active_by_side[side] + kChannelsPerChip - 1) / kChannelsPerChip;
  }

  estimate.active_channels += active_channels;
  estimate.packed_chips += packed_chips;
  estimate.payload_bits += static_cast<double>(active_channels) * kBitsPerHit * kSamplesPerEvent;
  estimate.packed_bits +=
      (static_cast<double>(packed_chips) * kOverheadBits +
       static_cast<double>(active_channels) * kBitsPerHit) *
      kSamplesPerEvent;
  estimate.per_channel_overhead_bits +=
      static_cast<double>(active_channels) * (kOverheadBits + kBitsPerHit) * kSamplesPerEvent;
}

void print_estimate(const char* name, const Estimate& estimate, std::uint64_t n_events) {
  const double total_time_sec = static_cast<double>(n_events) * kEventWindowSec;
  const double payload_rate_gbps = estimate.payload_bits / total_time_sec / 1e9;
  const double packed_rate_gbps = estimate.packed_bits / total_time_sec / 1e9;
  const double conservative_rate_gbps = estimate.per_channel_overhead_bits / total_time_sec / 1e9;

  std::cout << name << "\n"
            << "  mean active channels/event: "
            << static_cast<double>(estimate.active_channels) / static_cast<double>(n_events) << "\n"
            << "  mean packed chips/event:    "
            << static_cast<double>(estimate.packed_chips) / static_cast<double>(n_events) << "\n"
            << "  payload-only:              " << payload_rate_gbps << " Gb/s\n"
            << "  packed " << kChannelsPerChip << "-channel:         " << packed_rate_gbps << " Gb/s\n"
            << "  overhead per channel:      " << conservative_rate_gbps << " Gb/s\n"
            << "  payload underestimate:     "
            << 100.0 * (packed_rate_gbps - payload_rate_gbps) / packed_rate_gbps << "%\n"
            << "  per-channel overestimate:  "
            << 100.0 * (conservative_rate_gbps - packed_rate_gbps) / packed_rate_gbps << "%\n";
}

}  // namespace

void overhead_estimate(const char* input_dir = "data/test") {
  const auto files = find_root_files(input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << input_dir << "\n";
    return;
  }

  std::array<Estimate, 2> original_estimates;
  std::array<Estimate, 2> lfhcal_estimates;
  std::string collection_name;
  std::uint64_t n_events = 0;

  std::cerr << "Reading " << files.size() << " files.\n";
  for (std::size_t file_index = 0; file_index < files.size(); ++file_index) {
    std::cerr << "\r" << (file_index + 1) << "/" << files.size() << " files read." << std::flush;

    podio::ROOTReader reader;
    reader.openFile(files[file_index]);
    const std::size_t total_events = reader.getEntries("events");

    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (collection_name.empty()) collection_name = find_insert_collection(frame);
      if (collection_name.empty()) continue;
      const auto& collections = frame.getAvailableCollections();
      if (std::find(collections.begin(), collections.end(), collection_name) == collections.end()) continue;
      ++n_events;

      std::array<std::unordered_map<std::uint64_t, double>, 2> original_energies;
      std::array<std::unordered_map<std::uint64_t, double>, 2> original_no_ring_energies;
      std::array<std::unordered_map<VirtualChannel, double, VirtualChannelHash>, 2> lfhcal_energies;
      std::array<std::unordered_map<VirtualChannel, double, VirtualChannelHash>, 2> lfhcal_no_ring_energies;

      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(collection_name);
      for (const auto& hit : hits) {
        const std::uint64_t cell_id = static_cast<std::uint64_t>(hit.getCellID());
        const int side = decode_side(cell_id);
        const int layer = decode_layer(cell_id);
        if (side < 0 || side > 1 || layer < 1 || layer > 60) continue;

        const auto position = hit.getPosition();
        const double x_mm = static_cast<double>(position.x);
        const double y_mm = static_cast<double>(position.y);
        const double energy_gev = hit.getEnergy();
        original_energies[side][cell_id] += energy_gev;
        if (!inside_veto(x_mm, y_mm, kOriginalVetoRadiusMM)) {
          original_no_ring_energies[side][cell_id] += energy_gev;
        }

        const VirtualChannel channel{
            layer - 1,
            static_cast<int>(std::floor(x_mm / kVirtualCellSizeMM)),
            static_cast<int>(std::floor(y_mm / kVirtualCellSizeMM)),
        };
        lfhcal_energies[side][channel] += energy_gev;
        const double channel_x_mm = (static_cast<double>(channel.ix) + 0.5) * kVirtualCellSizeMM;
        const double channel_y_mm = (static_cast<double>(channel.iy) + 0.5) * kVirtualCellSizeMM;
        if (!inside_veto(channel_x_mm, channel_y_mm, kLFHCalVetoRadiusMM)) {
          lfhcal_no_ring_energies[side][channel] += energy_gev;
        }
      }

      std::array<std::uint64_t, 2> original_active{};
      std::array<std::uint64_t, 2> original_no_ring_active{};
      std::array<std::uint64_t, 2> lfhcal_active{};
      std::array<std::uint64_t, 2> lfhcal_no_ring_active{};
      for (int side = 0; side < 2; ++side) {
        for (const auto& [channel, energy] : original_energies[side]) {
          (void)channel;
          if (energy > kThresholdGeV) ++original_active[side];
        }
        for (const auto& [channel, energy] : original_no_ring_energies[side]) {
          (void)channel;
          if (energy > kThresholdGeV) ++original_no_ring_active[side];
        }
        for (const auto& [channel, energy] : lfhcal_energies[side]) {
          (void)channel;
          if (energy > kThresholdGeV) ++lfhcal_active[side];
        }
        for (const auto& [channel, energy] : lfhcal_no_ring_energies[side]) {
          (void)channel;
          if (energy > kThresholdGeV) ++lfhcal_no_ring_active[side];
        }
      }

      add_event(original_estimates[0], original_active);
      add_event(original_estimates[1], original_no_ring_active);
      add_event(lfhcal_estimates[0], lfhcal_active);
      add_event(lfhcal_estimates[1], lfhcal_no_ring_active);
    }
  }
  std::cerr << "\n";

  if (n_events == 0) {
    std::cerr << "No insert events found.\n";
    return;
  }

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "Events: " << n_events << "\n"
            << "Threshold: 0.5 MIP (" << kThresholdGeV << " GeV per single-layer channel)\n"
            << "Packing assumption: channels are packed independently on each side into groups of up to "
            << kChannelsPerChip << ".\n\n";
  print_estimate("Original insert, with inner ring", original_estimates[0], n_events);
  print_estimate("Original insert, without inner ring", original_estimates[1], n_events);
  print_estimate("LFHCal tiles, with inner ring", lfhcal_estimates[0], n_events);
  print_estimate("LFHCal tiles, without inner ring", lfhcal_estimates[1], n_events);
}
