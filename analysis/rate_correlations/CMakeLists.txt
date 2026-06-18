find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

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

install(TARGETS rate_scan rate_scan_r1000 rate_scan_origins_r1000 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
