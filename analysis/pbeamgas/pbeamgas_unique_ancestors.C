/*

root -l -b -q 'analysis/pbeamgas/pbeamgas_unique_ancestors.C("data_directory/path")'

*/

#include <TChain.h>
#include <TFile.h>
#include <TInterpreter.h>
#include <TROOT.h>
#include <TSystem.h>

#include <bit>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

// Find all files in a directory.
std::vector<std::string> find_root_files(const std::string& input_dir) {
  std::vector<std::string> files;
  for (const auto& entry : fs::recursive_directory_iterator(input_dir)) {
    if (!entry.is_regular_file()) continue;
    if (entry.path().extension() == ".root") files.push_back(entry.path().string());
  }
  std::sort(files.begin(), files.end());
  return files;
}

bool is_pbeamgas(int status) {
  return status >= 6000 && status < 7000;
}

// Create a hash for a MCParticle based on its PDG, generatorStatus, 4-vector, and vertex.
std::uint64_t particle_hash(const edm4hep::MCParticle& particle) {
  const auto momentum = particle.getMomentum();
  const auto vertex = particle.getVertex();
  std::uint64_t h = 0;
  const auto mix = [&](std::uint64_t value) {
    h ^= std::hash<std::uint64_t>{}(value) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  };
  mix(static_cast<std::uint64_t>(static_cast<std::uint32_t>(particle.getPDG())));
  mix(static_cast<std::uint64_t>(static_cast<std::uint32_t>(particle.getGeneratorStatus())));
  mix(std::bit_cast<std::uint64_t>(momentum.x));
  mix(std::bit_cast<std::uint64_t>(momentum.y));
  mix(std::bit_cast<std::uint64_t>(momentum.z));
  mix(std::bit_cast<std::uint64_t>(particle.getEnergy()));
  mix(std::bit_cast<std::uint64_t>(vertex.x));
  mix(std::bit_cast<std::uint64_t>(vertex.y));
  mix(std::bit_cast<std::uint64_t>(vertex.z));
  return h;
}

}  // namespace

void pbeamgas_unique_ancestors(const char* input_dir = "data/bkg_apr") {
  gInterpreter->Declare(R"CPP(
    #include <podio/Frame.h>
    #include <podio/ROOTReader.h>
    #include <edm4hep/MCParticle.h>
    #include <edm4hep/SimCalorimeterHitCollection.h>
  )CPP");

  const auto files = find_root_files(input_dir);
  if (files.empty()) {
    std::cerr << "No ROOT files found in " << input_dir << "\n";
    return;
  }

  std::unordered_map<std::uint64_t, std::uint64_t> particle_counts;
  std::uint64_t n_events = 0;
  std::uint64_t n_beamgas_contributions = 0;

  std::cerr << "Reading " << files.size() << " files.\n";

  // Loop through input files in directory.
  for (std::size_t file_index = 0; file_index < files.size(); ++file_index) {
    std::cerr << "\r" << (file_index + 1) << "/" << files.size() << " files read." << std::flush;

    podio::ROOTReader reader;
    reader.openFile(files[file_index]);
    const std::size_t total_events = reader.getEntries("events");

    // Loop all events in this file.
    for (std::size_t event_index = 0; event_index < total_events; ++event_index) {
      auto data = reader.readEntry("events", event_index);
      if (!data) continue;

      podio::Frame frame(std::move(data));
      const auto& collections = frame.getAvailableCollections();
      if (std::find(collections.begin(), collections.end(), "LFHCALHits") == collections.end()) continue;
      ++n_events;

      // Loop all hits in this event.
      const auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>("LFHCALHits");
      for (const auto& hit : hits) {

        // Loop all contributions in this hit.
        for (const auto& contribution : hit.getContributions()) {
          if (contribution.getEnergy() <= 0.0) continue;

          // Grab contrib particle
          auto particle = contribution.getParticle();

          // Keep only beam gas.
          if (!is_pbeamgas(particle.getGeneratorStatus())) continue;
          ++n_beamgas_contributions;

          // Walk up ancestry until no more parent particles.
          while (!particle.getParents().empty()) {
            particle = particle.getParents()[0];
          }

          // Create an identifying hash for the final ancestor particle.
          ++particle_counts[particle_hash(particle)];
        }
      }
    }
  }

  std::cout << "\n";
  std::cout << "Total events: " << n_events << "\n";
  std::cout << "Beam gas contributions: " << n_beamgas_contributions << "\n";
  std::cout << "Number of unique ancestors: " << particle_counts.size() << "\n";
}
