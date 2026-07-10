find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(insert_occupancy occupancy.cc)
target_compile_features(insert_occupancy PRIVATE cxx_std_17)
target_include_directories(insert_occupancy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(insert_occupancy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/insert_to_lfhcal.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/physical_layers/src/family.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/physical_layers/src/with_inner_ring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/physical_layers/src/no_inner_ring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/physical_layers/src/shared.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/physical_layers/src/xy.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/physical_layers/src/channel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/physical_layers/src/mip.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/family.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/with_inner_ring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/no_inner_ring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/shared.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/xy.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/channel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/chip.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/mip.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/radius.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lfhcal_segments/src/side.cc)
target_link_libraries(
  insert_occupancy
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

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

install(TARGETS insert_occupancy plot_hits insert_histogram_mip bkg_composition RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
