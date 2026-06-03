find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(EDM4EIC REQUIRED)

add_executable(occupancy occupancy.cc)
target_compile_features(occupancy PRIVATE cxx_std_17)
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
