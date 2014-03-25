cmake_minimum_required(VERSION 2.8)

project(creativity CXX)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

# No in-source building
include(MacroEnsureOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out-of-source build.  Create a build directory and run 'cmake ${CMAKE_SOURCE_DIR} [options]'.")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
#set(CMAKE_VERBOSE_MAKEFILE ON)

include(GNUInstallDirs)

# Use pkg-config to find the includes/libs we need
find_package(PkgConfig REQUIRED)
set(creativity_pkgconfig_GTKMM gtkmm-3.0)
set(creativity_pkgconfig_ERIS liberis)
set(creativity_pkgconfig ERIS GTKMM)
foreach(pkg ${creativity_pkgconfig})
    pkg_check_modules(${pkg} REQUIRED ${creativity_pkgconfig_${pkg}})
    include_directories(${${pkg}_INCLUDE_DIRS})
    link_directories(${${pkg}_LIBRARY_DIRS})
    add_definitions(${${pkg}_CFLAGS_OTHER})
endforeach()

#find_package(Boost REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

# Find library source/headers
file(GLOB_RECURSE cr_libsrc creativity/*.cpp)
file(GLOB_RECURSE cr_headers creativity/*.hpp)

# Set up compilation flags; put everything before CMAKE_CXX_FLAGS so that these can be overridden
set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -O3 -Wall -Wextra ${CMAKE_CXX_FLAGS}")
# Force on ERIS_DEBUG and -O0 for a Debug build
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -DERIS_DEBUG")

# Generate the config.hpp file (which needs to know where gui.glade was installed)
configure_file(${CMAKE_SOURCE_DIR}/creativity/config.hpp.in
    ${CMAKE_BINARY_DIR}/creativity/config.hpp)
include_directories(${CMAKE_BINARY_DIR})

add_library(creativity SHARED ${cr_libsrc})
foreach(pkg ${creativity_pkgconfig})
    target_link_libraries(creativity ${${pkg}_LIBRARIES})
endforeach()

add_executable(creativity-gui gui.cpp)
target_link_libraries(creativity-gui creativity)

install(FILES gui.glade creativity.svg DESTINATION ${CMAKE_INSTALL_DATADIR}/creativity/gui.glade)

#add_subdirectory(doc EXCLUDE_FROM_ALL)
