#pragma once

#include <podio/Frame.h>

#include <filesystem>
#include <string>
#include <vector>

namespace br {

bool has_collection(const podio::Frame& frame, const std::string& name);
std::vector<std::filesystem::path> find_root_files(const std::string& input_dir);

}  // namespace br
