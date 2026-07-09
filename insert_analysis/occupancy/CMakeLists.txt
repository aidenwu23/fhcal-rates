find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(insert_occupancy occupancy.cc)
target_compile_features(insert_occupancy PRIVATE cxx_std_17)
target_include_directories(insert_occupancy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_sources(insert_occupancy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/insert_to_lfhcal.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/shared.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/xy.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/channel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/chip.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mip.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/radius.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/side.cc)
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

install(TARGETS insert_occupancy RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
