cmake_minimum_required(VERSION 3.15)
project(NLOSingleInclusive CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

# Put the built binary directly in the top-level build directory rather
# than mirroring src/.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# --- GSL ---------------------------------------------------------------
find_package(GSL REQUIRED)

# --- LHAPDF --------------------------------------------------------------
# No CMake config is shipped by LHAPDF, so locate it manually.
# LHAPDF_ROOT can be set (-DLHAPDF_ROOT=/path/to/lhapdf) to override
# auto-detection, e.g. on machines where it isn't installed via Homebrew.
find_path(LHAPDF_INCLUDE_DIR LHAPDF/LHAPDF.h
    HINTS ${LHAPDF_ROOT}
    PATHS /opt/homebrew/opt/lhapdf /opt/homebrew/Cellar/lhapdf/*/  /usr/local/opt/lhapdf
    PATH_SUFFIXES include
)
find_library(LHAPDF_LIBRARY NAMES LHAPDF
    HINTS ${LHAPDF_ROOT}
    PATHS /opt/homebrew/opt/lhapdf /opt/homebrew/Cellar/lhapdf/*/  /usr/local/opt/lhapdf
    PATH_SUFFIXES lib
)

if(NOT LHAPDF_INCLUDE_DIR OR NOT LHAPDF_LIBRARY)
    message(FATAL_ERROR "LHAPDF not found. Set -DLHAPDF_ROOT=/path/to/lhapdf")
endif()

message(STATUS "Found LHAPDF include: ${LHAPDF_INCLUDE_DIR}")
message(STATUS "Found LHAPDF library: ${LHAPDF_LIBRARY}")

# The actual source list and target live in src/.
add_subdirectory(src)

target_include_directories(nlosingleinclusive PRIVATE ${LHAPDF_INCLUDE_DIR})
target_link_libraries(nlosingleinclusive PRIVATE GSL::gsl GSL::gslcblas ${LHAPDF_LIBRARY})
