cmake_minimum_required(VERSION 3.14)
project(DistributionTool VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Standard GNU install paths (bin, lib, share, …)
include(GNUInstallDirs)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Compiler flags.
# -march=native optimises for the build host CPU.
# For distributable/AppImage builds use: cmake … -DNATIVE_MARCH=OFF
option(NATIVE_MARCH "Optimise for the build host CPU (-march=native)" ON)
if(NATIVE_MARCH)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG -Wall")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall")
endif()
set(CMAKE_CXX_FLAGS_DEBUG    "-g -O0 -Wall -Wextra")

# ── Eigen3 (header-only) ────────────────────────────────────────────────────
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIR}")

# ── OpenMP (optional, speeds up neighbor-list & descriptor loops) ────────────
find_package(OpenMP)

# ── Library (shared between distool and tests) ────────────────────────────────
set(LIB_SOURCES
    src/LammpsDumpReader.cpp
    src/SphericalHarmonics.cpp
    src/RadialBasis.cpp
    src/SOAPDescriptor.cpp
    src/Statistics.cpp
    src/DefectClassifier.cpp
    src/PCA.cpp
)

add_library(distool_lib OBJECT ${LIB_SOURCES})
target_include_directories(distool_lib PUBLIC include)
target_link_libraries(distool_lib PUBLIC Eigen3::Eigen)

if(OpenMP_CXX_FOUND)
    target_link_libraries(distool_lib PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(distool_lib PUBLIC USE_OPENMP)
    message(STATUS "OpenMP support: enabled")
else()
    message(STATUS "OpenMP support: disabled")
endif()

# ── Main executable ───────────────────────────────────────────────────────────
add_executable(distool src/main.cpp)
target_link_libraries(distool PRIVATE distool_lib)

# ── Install rules (Linux) ─────────────────────────────────────────────────────
install(TARGETS distool
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Desktop integration (used by AppImage and system-wide installs)
install(FILES packaging/distool.desktop
        DESTINATION ${CMAKE_INSTALL_DATADIR}/applications)
install(FILES packaging/distool.png
        DESTINATION ${CMAKE_INSTALL_DATADIR}/icons/hicolor/256x256/apps)

# ── Tests ─────────────────────────────────────────────────────────────────────
enable_testing()
add_subdirectory(tests)

# ── Quick-build helper ────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "To build:")
message(STATUS "  cmake -B build -DCMAKE_BUILD_TYPE=Release")
message(STATUS "  cmake --build build -j$(nproc)")
message(STATUS "  ./build/distool --help")
