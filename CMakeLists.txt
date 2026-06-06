cmake_minimum_required(VERSION 3.16)

project(background-rates LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(scripts/data_prep)

add_subdirectory(analysis/hit-origins)

add_subdirectory(analysis/occupancy)

add_subdirectory(analysis/performance)
