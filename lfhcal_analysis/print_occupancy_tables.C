/*

root -l -b -q 'lfhcal_analysis/print_occupancy_tables.C("lfhcal_plots/occupancy.root")'

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

constexpr std::array<double, 3> kThresholds = {0.1, 0.3, 0.5};
constexpr double kChannelBitsPerEvent = 32.0 * 4.0;

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

std::string channel_rate(const TGraph& hit_graph, double threshold) {
  const double hit_rate_hz = graph_value(hit_graph, threshold);
  const double data_rate_gbps = hit_rate_hz * kChannelBitsPerEvent / 1.0e9;
  return format_significant(hit_rate_hz / 1.0e3, 3) + " kHz / " +
      format_significant(data_rate_gbps, 2) + " Gb/s";
}

std::string chip_rate(const TGraph& hit_graph, const TGraph& data_graph, double threshold) {
  return format_significant(graph_value(hit_graph, threshold) / 1.0e3, 3) + " kHz / " +
      format_significant(graph_value(data_graph, threshold), 2) + " Gb/s";
}

}  // namespace

void print_occupancy_tables(const char* input_file = "lfhcal_plots/occupancy.root",
                            const char* output_dir = "lfhcal_plots/tables") {
  // Open the consolidated occupancy output and locate every summary graph.
  TFile file(input_file, "READ");
  if (!file.IsOpen()) {
    std::cerr << "Failed to open " << input_file << "\n";
    return;
  }

  auto* channel_hit_p95 = get_graph(
      file, "occupancy/channel_hit_rate_vs_mip/graphs/g_channel_hit_rate_vs_mip_p95");
  auto* channel_hit_p99 = get_graph(
      file, "occupancy/channel_hit_rate_vs_mip/graphs/g_channel_hit_rate_vs_mip_p99");
  auto* chip_hit_p95 = get_graph(
      file, "occupancy/chip_hit_rate_vs_mip/graphs/g_chip_hit_rate_vs_mip_p95");
  auto* chip_hit_p99 = get_graph(
      file, "occupancy/chip_hit_rate_vs_mip/graphs/g_chip_hit_rate_vs_mip_p99");
  auto* chip_data_p95 = get_graph(
      file, "occupancy/chip_data_rate_vs_mip/graphs/g_chip_data_rate_vs_mip_p95");
  auto* chip_data_p99 = get_graph(
      file, "occupancy/chip_data_rate_vs_mip/graphs/g_chip_data_rate_vs_mip_p99");
  auto* total_lfhcal = get_graph(
      file, "occupancy/total_data_rate_vs_mip/graphs/g_lfhcal_rate");
  if (channel_hit_p95 == nullptr || channel_hit_p99 == nullptr ||
      chip_hit_p95 == nullptr || chip_hit_p99 == nullptr ||
      chip_data_p95 == nullptr || chip_data_p99 == nullptr ||
      total_lfhcal == nullptr) {
    return;
  }

  // Print the Markdown table and mirror its values into one CSV file.
  fs::create_directories(output_dir);
  std::ofstream csv(fs::path(output_dir) / "lfhcal.csv");
  csv << "MIP threshold,Channel p95 (hit/data),Channel p99 (hit/data),"
      << "Chip p95 (hit/data),Chip p99 (hit/data),Total LFHCAL rate\n";

  std::cout << "### Full LFHCAL\n\n"
            << "| MIP threshold | Channel p95 | Channel p99 | Chip p95 | Chip p99 | Total LFHCAL rate |\n"
            << "|---:|---:|---:|---:|---:|---:|\n";

  // Write one row for each representative threshold.
  for (double threshold : kThresholds) {
    const std::string channel_p95 = channel_rate(*channel_hit_p95, threshold);
    const std::string channel_p99 = channel_rate(*channel_hit_p99, threshold);
    const std::string chip_p95 = chip_rate(*chip_hit_p95, *chip_data_p95, threshold);
    const std::string chip_p99 = chip_rate(*chip_hit_p99, *chip_data_p99, threshold);
    const std::string total_data =
        format_significant(graph_value(*total_lfhcal, threshold), 3) + " Gb/s";

    std::cout << "| " << threshold << " | " << channel_p95 << " | " << channel_p99
              << " | " << chip_p95 << " | " << chip_p99 << " | " << total_data << " |\n";
    csv << threshold << "," << channel_p95 << "," << channel_p99 << ","
        << chip_p95 << "," << chip_p99 << "," << total_data << "\n";
  }
  std::cout << "\n";
}
