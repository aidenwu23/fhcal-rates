find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(lfhcal_occupancy occupancy/occupancy.cc)
target_compile_features(lfhcal_occupancy PRIVATE cxx_std_17)
target_include_directories(lfhcal_occupancy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy
  ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(lfhcal_occupancy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/contour.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/src/shared.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/src/family.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/src/xy.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/src/channel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/src/central.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/src/chip.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/src/mip.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/src/multiplicity.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/occupancy/src/radius.cc)
target_link_libraries(
  lfhcal_occupancy
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS lfhcal_occupancy RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(lfhcal_rate_summary rate_summary.cc)
target_compile_features(lfhcal_rate_summary PRIVATE cxx_std_17)
target_include_directories(lfhcal_rate_summary PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(lfhcal_rate_summary PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  lfhcal_rate_summary
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO)

add_executable(lfhcal_energy_spectrum energy_spectrum.cc)
target_compile_features(lfhcal_energy_spectrum PRIVATE cxx_std_17)
target_include_directories(lfhcal_energy_spectrum PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(lfhcal_energy_spectrum PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  lfhcal_energy_spectrum
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS lfhcal_rate_summary lfhcal_energy_spectrum RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(channel_count channel_count.cc)
target_compile_features(channel_count PRIVATE cxx_std_17)
target_include_directories(channel_count PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(channel_count PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  channel_count
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO)

install(TARGETS channel_count RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory(hit-origins)
add_subdirectory(muons)
