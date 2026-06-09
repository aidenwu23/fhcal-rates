find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(hit_origins origins.cc)
target_compile_features(hit_origins PRIVATE cxx_std_17)
target_include_directories(hit_origins PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_sources(hit_origins PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_channel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/eta.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/edep.cc)
target_link_libraries(
  hit_origins
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS hit_origins RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(eta_all eta_all.cc)
target_compile_features(eta_all PRIVATE cxx_std_17)
target_include_directories(eta_all PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(eta_all PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc)
target_link_libraries(
  eta_all
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS eta_all RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(multiplicity multiplicity.cc)
target_compile_features(multiplicity PRIVATE cxx_std_17)
target_include_directories(multiplicity PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(multiplicity PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc)
target_link_libraries(
  multiplicity
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad)

install(TARGETS multiplicity RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
