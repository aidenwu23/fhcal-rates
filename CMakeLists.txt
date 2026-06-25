cmake_minimum_required(VERSION 3.16)

project(background-rates LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

find_package(ROOT REQUIRED COMPONENTS Core RIO)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_subdirectory(scripts/data_prep)

add_subdirectory(analysis/hit-origins)

add_subdirectory(analysis/occupancy)

add_subdirectory(analysis/performance)

add_subdirectory(analysis/muons)

add_executable(channel_count analysis/channel_count.cc)
target_compile_features(channel_count PRIVATE cxx_std_17)
target_include_directories(channel_count PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(channel_count PRIVATE ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  channel_count
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO)

install(TARGETS channel_count RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
