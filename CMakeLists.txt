cmake_minimum_required(VERSION 2.8)

project(creativity CXX)

# Creativity package version
set(CREATIVITY_VERSION_MAJOR "0")
set(CREATIVITY_VERSION_MINOR "0")
set(CREATIVITY_VERSION_PATCH "1")
set(CREATIVITY_VERSION "${CREATIVITY_VERSION_MAJOR}.${CREATIVITY_VERSION_MINOR}.${CREATIVITY_VERSION_PATCH}")

set(creativity_author "Jason Rhinelander <jason@imaginary.ca>")
set(creativity_url "https://github.com/erisproject/creativity")
set(creativity_description "Agent-based model of creativity and piracy")

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
set(creativity_pkgconfig_PQXX libpqxx)

# pkg-config required packages for libcreativity.so:
set(creativity_pkgconfig ERIS)

option(DEBUG "Enable debugging output" OFF)
if (DEBUG)
    message(STATUS "Debugging enabled.")
    add_definitions("-DERIS_DEBUG")
else()
    message(STATUS "Debugging disabled.")
endif()

option(GUI "Enable gtkmm GUI interface" ON)
if (GUI)
    # required packages for libcreativitygui.so:
    set(creativity_gui_pkgconfig GTKMM)
else()
    message(STATUS "GUI disabled.")
    add_definitions("-DCREATIVITY_SKIP_GUI")
endif()

option(PGSQL "Enable PostgreSQL support" OFF)
if (PGSQL)
    list(APPEND creativity_pkgconfig PQXX)
else()
    message(STATUS "PostgreSQL support disabled.")
    add_definitions("-DCREATIVITY_SKIP_PGSQL")
endif()

option(STORAGE_THREAD "Enable threaded storage system" ON)
if (STORAGE_THREAD)
    message(STATUS "Threaded storage enabled.")
else()
    message(STATUS "Threaded storage disabled.")
    add_definitions("-DCREATIVITY_DISABLE_THREADED_STORAGE")
endif()

foreach(pkg ${creativity_pkgconfig} ${creativity_gui_pkgconfig})
    pkg_check_modules(${pkg} REQUIRED ${creativity_pkgconfig_${pkg}})
    include_directories(${${pkg}_INCLUDE_DIRS})
    link_directories(${${pkg}_LIBRARY_DIRS})
    add_definitions(${${pkg}_CFLAGS_OTHER})
endforeach()

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Find library source/headers
file(GLOB_RECURSE cr_libsrc creativity/*.cpp)
file(GLOB_RECURSE cr_gui_libsrc creativity/gui/*.cpp)
list(REMOVE_ITEM cr_libsrc ${cr_gui_libsrc})
if (NOT ${PGSQL})
    file(GLOB_RECURSE cr_psql creativity/state/PsqlStorage.cpp)
    list(REMOVE_ITEM cr_libsrc ${cr_psql})
endif()
# doc/CMakeLists.txt needs the list of headers:
file(GLOB_RECURSE cr_headers creativity/*.hpp)

# Set up compilation flags; put everything before CMAKE_CXX_FLAGS so that these can be overridden
set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -O3 -Wall -Wextra ${CMAKE_CXX_FLAGS}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # Increasing the inline threshold to 1000 makes a big difference for mt19937 RNG generation
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mllvm -inline-threshold=10000")
endif()

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

add_executable(creativity-cli cli.cpp)
target_link_libraries(creativity-cli creativity)

add_executable(creativity-info info.cpp)
target_link_libraries(creativity-info creativity)

add_executable(creativity-random random.cpp)
target_link_libraries(creativity-random ${ERIS_LIBRARIES})

if (GUI)
    add_library(creativitygui SHARED ${cr_gui_libsrc})
    foreach(pkg ${creativity_gui_pkgconfig})
        target_link_libraries(creativitygui ${${pkg}_LIBRARIES})
    endforeach()

    target_link_libraries(creativitygui creativity)

    add_executable(creativity-gui gui.cpp)
    target_link_libraries(creativity-gui creativitygui)

    # Copy the glade and icon to the install datadir
    install(FILES gui.glade creativity.svg DESTINATION ${CMAKE_INSTALL_DATADIR}/creativity/)

    # Also copy them into the build dir (so that running ./creativity-gui from build/ works)
    configure_file(gui.glade ${CMAKE_BINARY_DIR}/gui.glade COPYONLY)
    configure_file(creativity.svg ${CMAKE_BINARY_DIR}/creativity.svg COPYONLY)
endif()

add_subdirectory(doc)

add_subdirectory(scratch)
