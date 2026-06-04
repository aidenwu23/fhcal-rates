#include "utils.h"

#include <algorithm>

namespace br {

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

}  // namespace br
