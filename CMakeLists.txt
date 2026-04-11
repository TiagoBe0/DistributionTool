cmake_minimum_required(VERSION 3.14)
project(DistributionTool VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG    "-g -O0 -Wall -Wextra")

# ── Eigen3 (header-only) ────────────────────────────────────────────────────
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIR}")

# ── OpenMP (optional, speeds up neighbor-list & descriptor loops) ────────────
find_package(OpenMP)

# ── Sources ──────────────────────────────────────────────────────────────────
set(SOURCES
    src/LammpsDumpReader.cpp
    src/SphericalHarmonics.cpp
    src/RadialBasis.cpp
    src/SOAPDescriptor.cpp
    src/Statistics.cpp
    src/DefectClassifier.cpp
    src/PCA.cpp
    src/main.cpp
)

add_executable(distool ${SOURCES})
target_include_directories(distool PRIVATE include)
target_link_libraries(distool PRIVATE Eigen3::Eigen)

if(OpenMP_CXX_FOUND)
    target_link_libraries(distool PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(distool PRIVATE USE_OPENMP)
    message(STATUS "OpenMP support: enabled")
else()
    message(STATUS "OpenMP support: disabled")
endif()

install(TARGETS distool DESTINATION bin)

# ── Quick-build helper ────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "To build:")
message(STATUS "  cmake -B build -DCMAKE_BUILD_TYPE=Release")
message(STATUS "  cmake --build build -j$(nproc)")
message(STATUS "  ./build/distool --help")
