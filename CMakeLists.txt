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

add_subdirectory(insert_analysis)

add_subdirectory(lfhcal_analysis)
