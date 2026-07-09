find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_subdirectory(occupancy)
add_subdirectory(tile_channels)

add_executable(plot_hits plot_hits.cc)
target_compile_features(plot_hits PRIVATE cxx_std_17)
target_include_directories(plot_hits PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(plot_hits PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  plot_hits
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(bkg_composition bkg_composition.cc)
target_compile_features(bkg_composition PRIVATE cxx_std_17)
target_include_directories(bkg_composition PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(bkg_composition PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  bkg_composition
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(insert_histogram_mip histogram_mip.cc)
target_compile_features(insert_histogram_mip PRIVATE cxx_std_17)
target_include_directories(insert_histogram_mip PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(insert_histogram_mip PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  insert_histogram_mip
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS plot_hits insert_histogram_mip bkg_composition RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
