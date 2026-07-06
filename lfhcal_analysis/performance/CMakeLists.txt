find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(EDM4EIC REQUIRED)

add_executable(energy_res energy_res.cc)
target_compile_features(energy_res PRIVATE cxx_std_17)
target_include_directories(energy_res PRIVATE ${CMAKE_SOURCE_DIR}/scripts/include)
target_sources(energy_res PRIVATE ${CMAKE_SOURCE_DIR}/scripts/src/utils.cc)
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
  energy_res
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

install(TARGETS energy_res RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
