#pragma once

#include <podio/Frame.h>

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <vector>

namespace br {

bool has_collection(const podio::Frame& frame, const std::string& name);
std::vector<std::filesystem::path> find_root_files(const std::string& input_dir);
std::vector<double> log_edges(int bins, double low, double high);

class FileProgress {
 public:
  FileProgress(std::size_t total_files, std::ostream& os);
  void tick();

 private:
  std::size_t total_files_ = 0;
  std::size_t current_file_ = 0;
  std::ostream& os_;
};

}  // namespace br
