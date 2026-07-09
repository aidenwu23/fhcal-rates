find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(insert_tile_channels tile_channels.cc)
target_compile_features(insert_tile_channels PRIVATE cxx_std_17)
target_include_directories(insert_tile_channels PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(insert_tile_channels PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/insert_to_lfhcal.cc)
target_link_libraries(
  insert_tile_channels
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS insert_tile_channels RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
