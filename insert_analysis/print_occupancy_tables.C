/*

root -l -b -q 'insert_analysis/print_occupancy_tables.C("insert_plots/occupancy.root")'

*/

#include <TFile.h>
#include <TGraph.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

namespace fs = std::filesystem;

constexpr std::array<double, 3> kThresholds = {0.1, 0.5, 1.5};

struct Configuration {
  const char* directory;
  const char* title;
  bool has_chips;
};

constexpr std::array<Configuration, 3> kConfigurations = {{
    {"original_insert", "Original insert", true},
    {"lfhcal_tiles", "LFHCal-sized tiles", true},
    {"full_lfhcal", "Full LFHCal readout", true},
}};

struct RingConfiguration {
  const char* directory;
  const char* title;
};

constexpr std::array<RingConfiguration, 2> kRingConfigurations = {{
    {"with_inner_ring", "With inner ring"},
    {"no_inner_ring", "Without inner ring"},
}};

std::string format_significant(double value, int digits) {
  if (value == 0.0) return digits == 3 ? "0.00" : "0.0";

  const int exponent = static_cast<int>(std::floor(std::log10(std::abs(value))));
  const int decimals = std::max(0, digits - exponent - 1);
  const double scale = std::pow(10.0, static_cast<double>(digits - exponent - 1));
  const double rounded = std::round(value * scale) / scale;

  std::ostringstream text;
  text << std::fixed << std::setprecision(decimals) << rounded;
  return text.str();
}

TGraph* get_graph(TFile& file, const std::string& path) {
  auto* graph = dynamic_cast<TGraph*>(file.Get(path.c_str()));
  if (graph == nullptr) std::cerr << "Missing graph: " << path << "\n";
  return graph;
}

double graph_value(const TGraph& graph, double threshold) {
  for (int point = 0; point < graph.GetN(); ++point) {
    double x = 0.0;
    double y = 0.0;
    graph.GetPoint(point, x, y);
    if (std::abs(x - threshold) < 1e-9) return y;
  }
  return 0.0;
}

void print_table(TFile& file,
                 const Configuration& configuration,
                 const RingConfiguration& ring,
                 const fs::path& output_dir) {
  const std::string base = std::string(configuration.directory) + "/" + ring.directory;
  auto* hit_p95 = get_graph(file, base + "/channel_hit_rate_vs_mip/g_channel_hit_rate_vs_mip_p95");
  auto* hit_p99 = get_graph(file, base + "/channel_hit_rate_vs_mip/g_channel_hit_rate_vs_mip_p99");
  auto* data_p95 = get_graph(file, base + "/channel_data_rate_vs_mip/g_channel_data_rate_vs_mip_p95");
  auto* data_p99 = get_graph(file, base + "/channel_data_rate_vs_mip/g_channel_data_rate_vs_mip_p99");
  TGraph* chip_hit_p95 = nullptr;
  TGraph* chip_hit_p99 = nullptr;
  TGraph* chip_data_p95 = nullptr;
  TGraph* chip_data_p99 = nullptr;
  if (configuration.has_chips) {
    chip_hit_p95 = get_graph(file, base + "/chip_hit_rate_vs_mip/g_chip_hit_rate_vs_mip_p95");
    chip_hit_p99 = get_graph(file, base + "/chip_hit_rate_vs_mip/g_chip_hit_rate_vs_mip_p99");
    chip_data_p95 = get_graph(file, base + "/chip_data_rate_vs_mip/g_chip_data_rate_vs_mip_p95");
    chip_data_p99 = get_graph(file, base + "/chip_data_rate_vs_mip/g_chip_data_rate_vs_mip_p99");
  }
  auto* left = get_graph(file, base + "/side_data_rate_vs_mip/g_left_data_rate_vs_mip");
  auto* right = get_graph(file, base + "/side_data_rate_vs_mip/g_right_data_rate_vs_mip");
  if (hit_p95 == nullptr || hit_p99 == nullptr || data_p95 == nullptr ||
      data_p99 == nullptr || left == nullptr || right == nullptr) return;
  if (configuration.has_chips &&
      (chip_hit_p95 == nullptr || chip_hit_p99 == nullptr ||
       chip_data_p95 == nullptr || chip_data_p99 == nullptr)) return;

  std::cout << "### " << configuration.title << " — " << ring.title << "\n\n"
            << "| MIP threshold | Channel p95 | Channel p99 | Chip p95 | Chip p99 | Left | Right |\n"
            << "|---:|---:|---:|---:|---:|---:|---:|\n";

  const fs::path csv_path =
      output_dir / (std::string(configuration.directory) + "_" + ring.directory + ".csv");
  std::ofstream csv(csv_path);
  csv << "MIP threshold,Channel p95 (hit/data),Channel p99 (hit/data),"
      << "Chip p95 (hit/data),Chip p99 (hit/data),Left data rate,Right data rate\n";

  for (double threshold : kThresholds) {
    const std::string p95 =
        format_significant(graph_value(*hit_p95, threshold) / 1e3, 3) + " kHz / " +
        format_significant(graph_value(*data_p95, threshold), 2) + " Gb/s";
    const std::string p99 =
        format_significant(graph_value(*hit_p99, threshold) / 1e3, 3) + " kHz / " +
        format_significant(graph_value(*data_p99, threshold), 2) + " Gb/s";
    std::string chip_p95 = "—";
    std::string chip_p99 = "—";
    if (configuration.has_chips) {
      chip_p95 =
          format_significant(graph_value(*chip_hit_p95, threshold) / 1e3, 3) + " kHz / " +
          format_significant(graph_value(*chip_data_p95, threshold), 2) + " Gb/s";
      chip_p99 =
          format_significant(graph_value(*chip_hit_p99, threshold) / 1e3, 3) + " kHz / " +
          format_significant(graph_value(*chip_data_p99, threshold), 2) + " Gb/s";
    }
    const std::string left_rate = format_significant(graph_value(*left, threshold), 3) + " Gb/s";
    const std::string right_rate = format_significant(graph_value(*right, threshold), 3) + " Gb/s";

    std::cout << "| " << threshold << " | " << p95 << " | " << p99 << " | " << chip_p95 << " | " << chip_p99 << " | "
              << left_rate << " | " << right_rate << " |\n";
    csv << threshold << "," << p95 << "," << p99 << "," << chip_p95 << "," << chip_p99 << ","
        << left_rate << "," << right_rate << "\n";
  }
  std::cout << "\n";
}

}  // namespace

void print_occupancy_tables(const char* input_file = "insert_plots/occupancy.root",
                            const char* output_dir = "insert_plots/tables") {
  TFile file(input_file, "READ");
  if (!file.IsOpen()) {
    std::cerr << "Failed to open " << input_file << "\n";
    return;
  }

  fs::create_directories(output_dir);

  for (const auto& configuration : kConfigurations) {
    for (const auto& ring : kRingConfigurations) {
      print_table(file, configuration, ring, output_dir);
    }
  }
}
