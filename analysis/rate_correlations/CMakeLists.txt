find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(threshold_rate_overlays threshold_rate_overlays.cc)
target_compile_features(threshold_rate_overlays PRIVATE cxx_std_17)
target_include_directories(threshold_rate_overlays PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(threshold_rate_overlays PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  threshold_rate_overlays
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS threshold_rate_overlays RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
