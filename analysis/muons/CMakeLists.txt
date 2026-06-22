find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(channel_energy channel_energy.cc)
target_compile_features(channel_energy PRIVATE cxx_std_17)
target_include_directories(channel_energy PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(channel_energy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  channel_energy
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS channel_energy RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
