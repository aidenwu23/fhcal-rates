find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(chip_rate_vs_mip chip_rate_vs_mip.cc)
target_compile_features(chip_rate_vs_mip PRIVATE cxx_std_17)
target_include_directories(chip_rate_vs_mip PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(chip_rate_vs_mip PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  chip_rate_vs_mip
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(chip_rate_scan chip_rate_scan.cc)
target_compile_features(chip_rate_scan PRIVATE cxx_std_17)
target_include_directories(chip_rate_scan PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(chip_rate_scan PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  chip_rate_scan
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(fired_channels fired_channels.cc)
target_compile_features(fired_channels PRIVATE cxx_std_17)
target_include_directories(fired_channels PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(fired_channels PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  fired_channels
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(fired_channels_vs_rate fired_channels_vs_rate.cc)
target_compile_features(fired_channels_vs_rate PRIVATE cxx_std_17)
target_include_directories(fired_channels_vs_rate PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(fired_channels_vs_rate PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  fired_channels_vs_rate
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

add_executable(chip_channel_rate_vs_mip chip_channel_rate_vs_mip.cc)
target_compile_features(chip_channel_rate_vs_mip PRIVATE cxx_std_17)
target_include_directories(chip_channel_rate_vs_mip PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(chip_channel_rate_vs_mip PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc)
target_link_libraries(
  chip_channel_rate_vs_mip
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS chip_rate_vs_mip chip_rate_scan fired_channels fired_channels_vs_rate chip_channel_rate_vs_mip RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
