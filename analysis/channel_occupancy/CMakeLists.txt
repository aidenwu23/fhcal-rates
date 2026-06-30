find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(EDM4EIC REQUIRED)

add_executable(channel_occupancy channel_occupancy.cc)
target_compile_features(channel_occupancy PRIVATE cxx_std_17)
target_include_directories(channel_occupancy PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_sources(channel_occupancy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/smooth_hists.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/shared.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/reco.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/truth.cc)
set(_edm4eic_targets)
foreach(tgt EDM4EIC::edm4eic EDM4eic::edm4eic edm4eic::edm4eic)
  if(TARGET ${tgt})
    list(APPEND _edm4eic_targets ${tgt})
  endif()
endforeach()
if(NOT _edm4eic_targets)
  message(FATAL_ERROR "Could not determine EDM4eic link targets")
endif()

target_link_libraries(
  channel_occupancy
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ${_edm4eic_targets}
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS channel_occupancy RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(longitudinal_occupancy longitudinal_occupancy.cc)
target_compile_features(longitudinal_occupancy PRIVATE cxx_std_17)
target_include_directories(longitudinal_occupancy PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_sources(longitudinal_occupancy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/shared_longitudinal.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/reco_longitudinal.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/truth_longitudinal.cc)
target_link_libraries(
  longitudinal_occupancy
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ${_edm4eic_targets}
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS longitudinal_occupancy RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(channel_multiplicity channel_multiplicity.cc)
target_compile_features(channel_multiplicity PRIVATE cxx_std_17)
target_include_directories(channel_multiplicity PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(channel_multiplicity PRIVATE ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc)
target_link_libraries(
  channel_multiplicity
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ${_edm4eic_targets}
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS channel_multiplicity RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(rate_summary rate_summary.cc)
target_compile_features(rate_summary PRIVATE cxx_std_17)
target_include_directories(rate_summary PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(rate_summary PRIVATE ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  rate_summary
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO)

install(TARGETS rate_summary RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(energy_spectrum energy_spectrum.cc)
target_compile_features(energy_spectrum PRIVATE cxx_std_17)
target_include_directories(energy_spectrum PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(energy_spectrum PRIVATE ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  energy_spectrum
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS energy_spectrum RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(rate_scan rate_scan.cc)
target_compile_features(rate_scan PRIVATE cxx_std_17)
target_include_directories(rate_scan PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(rate_scan PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  rate_scan
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(rate_scan_r1000 rate_scan_r1000.cc)
target_compile_features(rate_scan_r1000 PRIVATE cxx_std_17)
target_include_directories(rate_scan_r1000 PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(rate_scan_r1000 PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  rate_scan_r1000
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(rate_scan_origins_r1000 rate_scan_origins_r1000.cc)
target_compile_features(rate_scan_origins_r1000 PRIVATE cxx_std_17)
target_include_directories(rate_scan_origins_r1000 PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(rate_scan_origins_r1000 PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc)
target_link_libraries(
  rate_scan_origins_r1000
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(rate_vs_mip rate_vs_mip.cc)
target_compile_features(rate_vs_mip PRIVATE cxx_std_17)
target_include_directories(rate_vs_mip PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(rate_vs_mip PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  rate_vs_mip
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(rate_vs_mip_r1000 rate_vs_mip_r1000.cc)
target_compile_features(rate_vs_mip_r1000 PRIVATE cxx_std_17)
target_include_directories(rate_vs_mip_r1000 PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(rate_vs_mip_r1000 PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  rate_vs_mip_r1000
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS rate_scan rate_scan_r1000 rate_scan_origins_r1000 rate_vs_mip rate_vs_mip_r1000 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
