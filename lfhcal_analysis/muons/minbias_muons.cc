/*

./build/minbias_muons -i data/bkg_july_minbias -o lfhcal_plots/muons/minbias_muons.root

*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1I.h>
#include <TH2D.h>

#include <podio/Frame.h>
#include <podio/ObjectID.h>
#include <podio/ROOTReader.h>

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include "decode_cell_id.h"
#include "utils.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kMCParticleCollection = "MCParticles";
constexpr const char* kHitCollection = "LFHCALHits";
constexpr int kNReadoutLayers = 8;
constexpr double kMinLFHCALEta = 1.2;
constexpr double kMaxLFHCALEta = 3.5;

// ----------------------------------------------------------------------------------
// CLI handling.
// ----------------------------------------------------------------------------------
struct Args {
  std::string input_path;
  std::string output_file;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT.root|DIR -o OUTPUT.root\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;

  // Loop arguments and read each path.
  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      args.input_path = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      args.output_file = argv[++i];
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }

  // Require both paths.
  if (args.input_path.empty() || args.output_file.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  return args;
}

// ----------------------------------------------------------------------------------
// Helpers.
// ----------------------------------------------------------------------------------
// Combine the collection and object indices into one particle identifier.
std::uint64_t object_key(const podio::ObjectID& id) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(id.collectionID)) << 32) |
         static_cast<std::uint32_t>(id.index);
}

// Build bin edges halfway between channel centers.
std::vector<double> make_edges(const std::set<double>& coordinates) {
  std::vector<double> values(coordinates.begin(), coordinates.end());
  if (values.empty()) return {};
  if (values.size() == 1) return {values.front() - 0.5, values.front() + 0.5};

  // Place interior edges between adjacent centers.
  std::vector<double> edges(values.size() + 1);
  edges.front() = values.front() - 0.5 * (values[1] - values[0]);
  for (std::size_t i = 1; i < values.size(); ++i) {
    edges[i] = 0.5 * (values[i - 1] + values[i]);
  }
  edges.back() = values.back() + 0.5 * (values.back() - values[values.size() - 2]);
  return edges;
}

}  // namespace

// ----------------------------------------------------------------------------------
// Main.
// ----------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  // Keep histograms independent of the output file and include overflow statistics.
  TH1::AddDirectory(false);
  TH1::StatOverflows(true);

  // Find input ROOT files.
  const auto args = parse_args(argc, argv);
  const fs::path input_path(args.input_path);
  std::vector<fs::path> files;
  if (fs::is_regular_file(input_path) && input_path.extension() == ".root") {
    files.push_back(input_path);
  } else {
    files = rates::find_root_files(args.input_path);
  }
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << args.input_path << "\n";
    return 1;
  }

  // Create the output directory.
  const fs::path output_path(args.output_file);
  if (output_path.has_parent_path()) fs::create_directories(output_path.parent_path());

  // Prepare event-level histograms.
  TH1I h_muons_event("h_muons_event", "Stored MC muons;muons/event;events", 200, -0.5, 199.5);
  TH1I h_lfhcal_direction_muons_event("h_lfhcal_direction_muons_event", "Stored MC muons with 1.2 < #eta < 3.5;muons/event;events", 200, -0.5, 199.5);
  TH1I h_lfhcal_muons_event("h_lfhcal_muons_event", "Muons depositing energy in LFHCAL;muons/event;events", 200, -0.5, 199.5);

  // Prepare sample totals and per-channel map inputs.
  const rates::LFHCALCellIDDecoder decoder;
  std::uint64_t events = 0;
  std::uint64_t muons = 0;
  std::uint64_t lfhcal_direction_muons = 0;
  std::uint64_t lfhcal_muons = 0;
  std::array<std::unordered_map<rates::LFHCALChannelID,
                                rates::LFHCALCellPosition,
                                rates::LFHCALChannelIDHash>,
             kNReadoutLayers> channel_positions;
  std::array<std::unordered_map<rates::LFHCALChannelID,
                                std::uint64_t,
                                rates::LFHCALChannelIDHash>,
             kNReadoutLayers> channel_muon_counts;
  rates::FileProgress progress(files.size(), std::cerr);

  // Loop files.
  for (const auto& path : files) {
    progress.tick();

    // Open the file and read its event count.
    podio::ROOTReader reader;
    reader.openFile(path.string());
    const std::size_t entries = reader.getEntries("events");

    // Loop events.
    for (std::size_t event_index = 0; event_index < entries; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      if (!rates::has_collection(frame, kMCParticleCollection) ||
          !rates::has_collection(frame, kHitCollection)) {
        continue;
      }
      ++events;

      // Collect all stored muons and those inside the LFHCal acceptance.
      std::unordered_set<std::uint64_t> event_muons;
      std::uint64_t event_lfhcal_direction_muons = 0;
      const auto& particles = frame.get<edm4hep::MCParticleCollection>(kMCParticleCollection);

      // Loop MC particles and filter by pseudorapidity.
      for (const auto& particle : particles) {
        if (particle.getPDG() != 13 && particle.getPDG() != -13) continue;
        event_muons.insert(object_key(particle.getObjectID()));
        const auto momentum = particle.getMomentum();
        const double eta = std::asinh(momentum.z / std::hypot(momentum.x, momentum.y));

        if (eta > kMinLFHCALEta && eta < kMaxLFHCALEta) ++event_lfhcal_direction_muons;
      }

      // Collect LFHCal muons once per event and once per channel.
      std::unordered_set<std::uint64_t> event_lfhcal_muons;
      std::array<std::unordered_map<rates::LFHCALChannelID,
                                    std::unordered_set<std::uint64_t>,
                                    rates::LFHCALChannelIDHash>,
                 kNReadoutLayers> layer_channel_muons;
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>(kHitCollection);

      // Loop hits.
      for (const auto& hit : hits) {
        // Decode active hits and record their channel centers.
        const auto cell_id = static_cast<std::uint64_t>(hit.getCellID());
        if (decoder.is_passive(cell_id)) continue;

        const auto channel = decoder.channel(cell_id);
        if (channel.rlayerz < 0 || channel.rlayerz >= kNReadoutLayers) continue;
        channel_positions[channel.rlayerz][channel] = decoder.position(cell_id);

        // Loop contributions.
        for (const auto& contribution : hit.getContributions()) {
          
          // Keep positive deposits linked to stored muons.
          if (contribution.getEnergy() <= 0.0) continue;
          const auto particle = contribution.getParticle();
          if (!particle.isAvailable()) continue;
          if (particle.getPDG() != 13 && particle.getPDG() != -13) continue;

          const auto id = object_key(particle.getObjectID());
          if (event_muons.count(id) == 0) continue;

          // Record the event-level muon and its channel crossing.
          event_lfhcal_muons.insert(id);
          layer_channel_muons[channel.rlayerz][channel].insert(id);
        }
      }

      // Count each muon once per channel.
      for (int layer = 0; layer < kNReadoutLayers; ++layer) {
        for (const auto& [channel, channel_muons] : layer_channel_muons[layer]) {
          channel_muon_counts[layer][channel] += channel_muons.size();
        }
      }

      // Fill event distributions and sample totals.
      h_muons_event.Fill(event_muons.size());
      h_lfhcal_direction_muons_event.Fill(event_lfhcal_direction_muons);
      h_lfhcal_muons_event.Fill(event_lfhcal_muons.size());
      muons += event_muons.size();
      lfhcal_direction_muons += event_lfhcal_direction_muons;
      lfhcal_muons += event_lfhcal_muons.size();
    }
  }
  std::cerr << "\n";

  // Require at least one processed event.
  if (events == 0) {
    std::cerr << "No readable events found\n";
    return 1;
  }

  // Build common X-Y axes from observed channel centers.
  std::set<double> x_coordinates;
  std::set<double> y_coordinates;
  // Loop layers and channels.
  for (const auto& layer_positions : channel_positions) {
    for (const auto& [channel, position] : layer_positions) {
      (void)channel;
      x_coordinates.insert(position.x_mm);
      y_coordinates.insert(position.y_mm);
    }
  }
  const auto x_edges = make_edges(x_coordinates);
  const auto y_edges = make_edges(y_coordinates);
  // Require spatial coordinates for the maps.
  if (x_edges.empty() || y_edges.empty()) {
    std::cerr << "No active LFHCAL channels found\n";
    return 1;
  }

  // Prepare one X-Y map per readout layer.
  std::array<TH2D*, kNReadoutLayers> h_muon_xy{};
  // Loop layers.
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    // Use the common channel-aligned axes.
    h_muon_xy[layer] = new TH2D(
        ("h_muon_xy_layer" + std::to_string(layer)).c_str(),
        ("LFHCAL muons, readout layer " + std::to_string(layer) + ";x [mm];y [mm];muons").c_str(),
        static_cast<int>(x_edges.size()) - 1, x_edges.data(),
        static_cast<int>(y_edges.size()) - 1, y_edges.data());
    // Loop channels reached by muons.
    for (const auto& [channel, count] : channel_muon_counts[layer]) {
      const auto position = channel_positions[layer].at(channel);
      h_muon_xy[layer]->Fill(position.x_mm, position.y_mm, count);
    }
  }

  // Write ROOT output.
  TFile output(args.output_file.c_str(), "RECREATE");
  // Stop if the ROOT file cannot be created.
  if (!output.IsOpen()) {
    std::cerr << "Failed to open output " << args.output_file << "\n";
    return 1;
  }

  // Write event-level histograms.
  h_muons_event.Write();
  h_lfhcal_direction_muons_event.Write();
  h_lfhcal_muons_event.Write();
  // Loop layers and draw X-Y maps.
  for (int layer = 0; layer < kNReadoutLayers; ++layer) {
    h_muon_xy[layer]->Write();

    // Draw a square map with space for the color scale.
    TCanvas canvas(("c_muon_xy_layer" + std::to_string(layer)).c_str(),
                   h_muon_xy[layer]->GetTitle(), 1000, 900);
    canvas.SetLeftMargin(0.12);
    canvas.SetRightMargin(0.16);
    canvas.SetTopMargin(0.10);
    canvas.SetBottomMargin(0.10);
    canvas.SetLogz();
    h_muon_xy[layer]->SetStats(false);
    h_muon_xy[layer]->Draw("COLZ");
    canvas.Write();
  }
  output.Close();

  // Write sample means.
  fs::path csv_path = output_path;
  csv_path.replace_extension(".csv");
  std::ofstream csv(csv_path);
  // Stop if the CSV file cannot be created.
  if (!csv) {
    std::cerr << "Failed to open output " << csv_path << "\n";
    return 1;
  }
  // Write one row of per-event means.
  csv << "events,muons_per_event,lfhcal_direction_muons_per_event,lfhcal_muons_per_event\n";
  csv << std::setprecision(10) << events << ','
      << static_cast<double>(muons) / events << ','
      << static_cast<double>(lfhcal_direction_muons) / events << ','
      << static_cast<double>(lfhcal_muons) / events << '\n';

  return 0;
}
