find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(EDM4EIC REQUIRED)

add_executable(occupancy occupancy.cc)
target_compile_features(occupancy PRIVATE cxx_std_17)
target_include_directories(occupancy PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_sources(occupancy PRIVATE
  ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/decode_cell_id.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/smooth_hists.cc
  ${CMAKE_SOURCE_DIR}/scripts/src/classify_hit.cc
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
  occupancy
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

install(TARGETS occupancy RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

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
