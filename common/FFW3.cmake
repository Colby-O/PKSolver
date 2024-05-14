message("-- PARAMETER_EXTRACTION - Loading FFW3...")

find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)

include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)