find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(data_rate_vs_mip data_rate_vs_mip.cc)
target_compile_features(data_rate_vs_mip PRIVATE cxx_std_17)
target_include_directories(data_rate_vs_mip PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(data_rate_vs_mip PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  data_rate_vs_mip
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(data_rate_scan data_rate_scan.cc)
target_compile_features(data_rate_scan PRIVATE cxx_std_17)
target_include_directories(data_rate_scan PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(data_rate_scan PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  data_rate_scan
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(data_rate_vs_radius data_rate_vs_radius.cc)
target_compile_features(data_rate_vs_radius PRIVATE cxx_std_17)
target_include_directories(data_rate_vs_radius PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(data_rate_vs_radius PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  data_rate_vs_radius
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(lfhcal_rate lfhcal_rate.cc)
target_compile_features(lfhcal_rate PRIVATE cxx_std_17)
target_include_directories(lfhcal_rate PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(lfhcal_rate PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  lfhcal_rate
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS data_rate_vs_mip data_rate_scan data_rate_vs_radius lfhcal_rate RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
