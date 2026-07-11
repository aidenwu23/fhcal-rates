find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(histogram_mip histogram_mip.cc)
target_compile_features(histogram_mip PRIVATE cxx_std_17)
target_include_directories(histogram_mip PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(histogram_mip PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  histogram_mip
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS histogram_mip RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
