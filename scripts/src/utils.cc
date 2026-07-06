#include "utils.h"

#include <algorithm>
#include <cmath>
#include <ostream>

namespace rates {

namespace {

constexpr double MIP_1 = 3.5e-3;
constexpr double MIP_2 = 7.0e-3;

}  // namespace

bool has_collection(const podio::Frame& frame, const std::string& name) {
  for (const auto& available : frame.getAvailableCollections()) {
    if (available == name) return true;
  }
  return false;
}

std::vector<std::filesystem::path> find_root_files(const std::string& input_dir) {
  std::vector<std::filesystem::path> files;
  for (const auto& entry : std::filesystem::directory_iterator(input_dir)) {
    if (entry.is_regular_file() && entry.path().extension() == ".root") {
      files.push_back(entry.path());
    }
  }
  std::sort(files.begin(), files.end());
  return files;
}

std::vector<double> log_edges(int bins, double low, double high) {
  std::vector<double> edges(bins + 1);
  const double log_low = std::log10(low);
  const double log_high = std::log10(high);
  for (int i = 0; i <= bins; ++i) {
    const double frac = static_cast<double>(i) / bins;
    edges[i] = std::pow(10.0, log_low + frac * (log_high - log_low));
  }
  return edges;
}

double mip_energy_gev(int layer) {
  return layer < 2 ? MIP_1 : MIP_2;
}

FileProgress::FileProgress(std::size_t total_files, std::ostream& os)
    : total_files_(total_files), os_(os) {}

void FileProgress::tick() {
  ++current_file_;
  os_ << "\r" << current_file_ << "/" << total_files_ << " files read." << std::flush;
}

}  // namespace rates
