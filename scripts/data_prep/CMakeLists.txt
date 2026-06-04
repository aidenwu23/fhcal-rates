find_package(ROOT REQUIRED COMPONENTS Core RIO)
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)

add_executable(filter filter.cc)
target_compile_features(filter PRIVATE cxx_std_17)
target_link_libraries(
  filter
  PRIVATE
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    ROOT::Core
    ROOT::RIO)

install(TARGETS filter
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
