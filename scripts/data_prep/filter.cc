/*

./build/filter \
  -i data/raw/epic/FULL/26.04.1/epic_craterlake/Bkg_Exact1S_2us/GoldCt/10um/DIS/NC/10x275/minQ2=1/pythia8NCDIS_10x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.0000.edm4hep.root \
  -o data/filtered/lfhcal_0000.root

*/

#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <podio/ROOTWriter.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct Args {
  std::string input;
  std::string output;
};

void usage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " -i INPUT.root -o OUTPUT.root\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string_view v(argv[i]);
    if ((v == "-i" || v == "--input") && i + 1 < argc) {
      args.input = argv[++i];
    } else if ((v == "-o" || v == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else {
      usage(argv[0]);
      std::exit(1);
    }
  }
  if (args.input.empty() || args.output.empty()) {
    usage(argv[0]);
    std::exit(1);
  }
  return args;
}

bool keepCollection(const std::string& name) {
  return name == "EventHeader" || name == "MCParticles" || name.find("LFHCAL") != std::string::npos ||
         name.find("HcalEndcapPInsert") != std::string::npos;
}

std::vector<std::string> keptCollections(const podio::Frame& frame) {
  std::vector<std::string> kept;
  for (const auto& name : frame.getAvailableCollections()) {
    if (keepCollection(name)) {
      kept.emplace_back(name);
    }
  }
  return kept;
}

} // namespace

int main(int argc, char* argv[]) {
  const Args args = parse_args(argc, argv);

  podio::ROOTReader reader;
  reader.openFile(args.input);

  const std::string category = "events";
  const auto entries = reader.getEntries(category);

  podio::ROOTWriter writer(args.output);

  std::size_t written = 0;
  for (std::size_t i = 0; i < entries; ++i) {
    auto data = reader.readEntry(category, i);
    podio::Frame frame(std::move(data));
    const auto kept = keptCollections(frame);
    writer.writeFrame(frame, category, kept);
    ++written;
  }

  writer.finish();
  std::cout << "wrote " << written << " events to " << args.output << '\n';
  return 0;
}
