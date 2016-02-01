cmake_minimum_required(VERSION 2.8)

project(creativity CXX)

# Creativity package version
set(CREATIVITY_VERSION_MAJOR "1")
set(CREATIVITY_VERSION_MINOR "2")
set(CREATIVITY_VERSION_PATCH "0")
set(CREATIVITY_VERSION "${CREATIVITY_VERSION_MAJOR}.${CREATIVITY_VERSION_MINOR}.${CREATIVITY_VERSION_PATCH}")

set(creativity_author "Jason Rhinelander <jason@imaginary.ca>")
set(creativity_url "https://imaginary.ca/eris/creativity")
set(creativity_description "Agent-based model of creativity and piracy")

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

set(ERIS_VERSION_REQUIRED "0.3.5")

# library versions (CURRENT.REVISION.AGE), which is totally separate
# from the above.  See http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
# Basic rules for updating these when releasing a new eris version:
# - If the new release has an API change:
#   - CURRENT gets updated
#   - REVISION gets set to 0
#   - If the API change was purely an addition, increment AGE
#   - Otherwise (i.e. an API change or removal), reset AGE to 0
# - Otherwise, if the release has no API change but has code changes:
#   - REVISION gets incremented
#   - (CURRENT and AGE stay the same)
# If there is no code change (e.g. the release is just a documentation update)
# then none of these change.
#
# (So something like 3.7.1 indicates the 8th revision of the libwhatever-3
# interface, and that code that links against liberis-2.* can link against this
# version, but code that links against liberis-1.* cannot.
# Main library so:
set(LIBCR_CUR 0)
set(LIBCR_REV 0)
set(LIBCR_AGE 0)
# Data library so:
set(LIBCRD_CUR 0)
set(LIBCRD_REV 0)
set(LIBCRD_AGE 0)
# GUI library so:
set(LIBCRG_CUR 0)
set(LIBCRG_REV 0)
set(LIBCRG_AGE 0)

# No in-source building
include(MacroEnsureOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out-of-source build.  Create a build directory and run 'cmake ${CMAKE_SOURCE_DIR} [options]'.")

#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
#set(CMAKE_VERBOSE_MAKEFILE ON)

include(GNUInstallDirs)

include(CheckIncludeFileCXX)

# Use pkg-config to find the includes/libs we need
find_package(PkgConfig REQUIRED)
set(creativity_pkgconfig_GTKMM gtkmm-3.0)
set(creativity_pkgconfig_ERIS liberis)
set(creativity_pkgconfig_LIBPNG libpng)
set(creativity_pkgconfig_CAIRO cairomm-1.0 cairomm-svg-1.0 cairomm-pdf-1.0)
set(creativity_pkgconfig_PANGO pangomm-1.4)

# pkg-config required packages for libcreativity.so:
set(creativity_pkgconfig ERIS)

# pkg-config required packages for libcreativity-data.so:
set(creativity_data_pkgconfig)

option(DEBUG "Enable debugging output" OFF)
if (DEBUG)
    message(STATUS "Debugging enabled (-DDEBUG=off to disable)")
    add_definitions("-DERIS_DEBUG")
else()
    message(STATUS "Debugging disabled (-DDEBUG=on to enable)")
endif()

option(GUI "Enable gtkmm GUI interface" ON)
if (GUI)
    # required packages for libcreativity-gui.so:
    set(creativity_gui_pkgconfig GTKMM)
    message(STATUS "GUI enabled (-DGUI=off to disable)")
else()
    message(STATUS "GUI disabled (-DGUI=on to enable)")
    add_definitions("-DCREATIVITY_SKIP_GUI")
endif()

option(DATA "Enable data parsing and results library and binaries" ON)
option(DATA_GRAPHING "Enable data graphing capabilites (requires DATA=ON)" ON)
if (DATA)
    message(STATUS "Data libraries and binaries enabled (-DDATA=off to disable)")
    if (DATA_GRAPHING)
        message(STATUS "    ... with graphing enabled (-DDATA_GRAPHING=OFF to disable)")
        set(creativity_data_pkgconfig ${creativity_data_pkgconfig} CAIRO PANGO LIBPNG)
    else()
        message(STATUS "    ... without graphing enabled (-DDATA_GRAPHING=ON to enable)")
    endif()
else()
    message(STATUS "Data libraries and binaries disabled (-DDATA=on to enable)")
endif()

option(STORAGE_THREAD "Enable threaded storage system" ON)
if (STORAGE_THREAD)
    message(STATUS "Threaded storage enabled (-DSTORAGE_THREAD=off to disable)")
else()
    message(STATUS "Threaded storage disabled (-DSTORAGE_THREAD=on to enable)")
    add_definitions("-DCREATIVITY_DISABLE_THREADED_STORAGE")
endif()

# Look for include-what-you-use, if requested
option(IWYU "Enable include-what-you-use header debugging support" OFF)
if (IWYU)
    cmake_minimum_required(VERSION 3.3)
    find_program(IWYU_PATH NAMES "iwyu" "include-what-you-use")
    if (NOT IWYU_PATH)
        message(FATAL_ERROR "Could not find executable for include-what-you-use support")
    endif()
    set(IWYU_ARGS ${IWYU_PATH} "-Xiwyu" "--mapping_file=${CMAKE_SOURCE_DIR}/creativity.imp")
    message(STATUS "include-what-you-use header debugging support enabled (-DIWYU=0 to disable)")
else()
    message(STATUS "include-what-you-use header debugging support disabled (-DIWYU=1 to enable)")
endif()

option(API_DOCS "Enable creativity API documentation generation (requires doxygen)" ON)
if (API_DOCS)
    message(STATUS "API documentation generation enabled (-DAPI_DOCS=off to disable)")
else()
    message(STATUS "API documentation generation disabled (-DAPI_DOCS=on to enable)")
endif()

foreach(pkg ${creativity_pkgconfig} ${creativity_gui_pkgconfig} ${creativity_data_pkgconfig})
    pkg_check_modules(${pkg} REQUIRED ${creativity_pkgconfig_${pkg}})
    include_directories(${${pkg}_INCLUDE_DIRS})
    link_directories(${${pkg}_LIBRARY_DIRS})
    add_definitions(${${pkg}_CFLAGS_OTHER})
endforeach()

if (DATA AND DATA_GRAPHING)
    CHECK_INCLUDE_FILE_CXX("png++/png.hpp" PNGPP_FOUND)
    if (NOT PNGPP_FOUND)
        unset(PNGPP_FOUND CACHE)
        message(FATAL_ERROR "libpng++ not found")
    endif()
endif()

if(ERIS_VERSION VERSION_LESS "${ERIS_VERSION_REQUIRED}")
    message(FATAL_ERROR "Eris v${ERIS_VERSION_REQUIRED} or above is required, but only v${ERIS_VERSION} was found")
endif()

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost REQUIRED COMPONENTS filesystem program_options system)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# Find library source/headers
file(GLOB_RECURSE cr_libsrc creativity/*.cpp)
file(GLOB_RECURSE cr_data_libsrc creativity/data/*.cpp creativity/cmdargs/Data.cpp creativity/cmdargs/Series*.cpp creativity/cmdargs/Results.cpp)
file(GLOB_RECURSE cr_gui_libsrc creativity/gui/*.cpp creativity/cmdargs/GUI.cpp)
list(REMOVE_ITEM cr_libsrc ${cr_gui_libsrc} ${cr_data_libsrc})
# doc/CMakeLists.txt needs the list of headers:
file(GLOB_RECURSE cr_headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} creativity/*.hpp creativity/*.hpp.in)

if (NOT GUI)
    # When not building the GUI, also omit the GUI headers
    file(GLOB_RECURSE cr_gui_headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} creativity/gui/*.hpp creativity/cmdargs/GUI.hpp)
    list(REMOVE_ITEM cr_headers ${cr_gui_headers})
endif()

if (NOT DATA)
    # When not building the data library/binaries, also omit the data headers
    file(GLOB_RECURSE cr_data_headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} creativity/data/*.hpp creativity/cmdargs/Data.hpp creativity/cmdargs/Series*.hpp creativity/cmdargs/Results.hpp)
    list(REMOVE_ITEM cr_headers ${cr_data_headers})
elseif(NOT DATA_GRAPHING)
    # Remove the data graphing headers
    file(GLOB_RECURSE cr_data_graph_headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} creativity/data/graph/*.hpp creativity/cmdargs/SeriesGraphs.hpp)
    list(REMOVE_ITEM cr_headers ${cr_data_graph_headers})
    # Remove the source files from the data src list:
    file(GLOB_RECURSE cr_data_graph_sources ${CMAKE_CURRENT_SOURCE_DIR} creativity/data/graph/*.cpp creativity/cmdargs/SeriesGraphs.cpp)
    list(REMOVE_ITEM cr_data_libsrc ${cr_data_graph_sources})
endif()

foreach(h ${cr_headers})
    if(h MATCHES "\\.in$")
        string(REGEX REPLACE "\\.in$" "" h_processed "${h}")
        configure_file("${h}" "${h_processed}" @ONLY)
        list(APPEND creativity_doc_include "${CMAKE_CURRENT_BINARY_DIR}/${h_processed}")
    else()
        list(APPEND creativity_doc_include "${CMAKE_CURRENT_SOURCE_DIR}/${h}")
    endif()
endforeach()
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Set up compilation flags; put everything before CMAKE_CXX_FLAGS so that these can be overridden
set(CMAKE_CXX_FLAGS "-std=c++11 -O2 -Wall -Wextra ${CMAKE_CXX_FLAGS}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # Increasing the inline threshold to 1000 makes a big difference for mt19937 RNG generation
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mllvm -inline-threshold=1000")
endif()

# Add -g for a Debug build
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

add_library(libcreativity SHARED ${cr_libsrc})
set_target_properties(libcreativity PROPERTIES
    OUTPUT_NAME creativity
    VERSION "${LIBCR_CUR}.${LIBCR_REV}.${LIBCR_AGE}"
    SOVERSION "${LIBCR_CUR}"
)
foreach(pkg ${creativity_pkgconfig})
    target_link_libraries(libcreativity ${${pkg}_LIBRARIES})
endforeach()
target_link_libraries(libcreativity ${Boost_LIBRARIES})
list(APPEND creativity_install_targets libcreativity)

add_executable(creativity-cli cli.cpp)
target_link_libraries(creativity-cli libcreativity)
list(APPEND creativity_install_targets creativity-cli)

add_executable(creativity-random random.cpp)
target_link_libraries(creativity-random ${ERIS_LIBRARIES})
list(APPEND creativity_install_targets creativity-random)

add_executable(creativity-info info.cpp)
target_link_libraries(creativity-info libcreativity)
list(APPEND creativity_install_targets creativity-info)

if (DATA)
    add_library(libcreativity-data SHARED ${cr_data_libsrc})
    set_target_properties(libcreativity-data PROPERTIES
        OUTPUT_NAME creativity-data
        VERSION "${LIBCRD_CUR}.${LIBCRD_REV}.${LIBCRD_AGE}"
        SOVERSION "${LIBCRD_CUR}"
    )
    foreach(pkg ${creativity_data_pkgconfig})
        target_link_libraries(libcreativity-data ${${pkg}_LIBRARIES})
    endforeach()

    target_link_libraries(libcreativity-data libcreativity)
    list(APPEND creativity_install_targets libcreativity-data)

    add_executable(creativity-data data.cpp)
    target_link_libraries(creativity-data libcreativity-data)
    list(APPEND creativity_install_targets creativity-data)

    add_executable(creativity-results results.cpp)
    target_link_libraries(creativity-results libcreativity-data)
    list(APPEND creativity_install_targets creativity-results)

    add_executable(creativity-series series.cpp)
    target_link_libraries(creativity-series libcreativity-data)
    list(APPEND creativity_install_targets creativity-series)

    add_executable(creativity-series-quantiles series-quantiles.cpp)
    target_link_libraries(creativity-series-quantiles libcreativity-data)
    list(APPEND creativity_install_targets creativity-series-quantiles)

    if (DATA_GRAPHING)
        add_executable(creativity-series-graphs series-graphs.cpp)
        target_link_libraries(creativity-series-graphs libcreativity-data)
        list(APPEND creativity_install_targets creativity-series-graphs)
    endif()
endif()


if (GUI)
    add_library(libcreativity-gui SHARED ${cr_gui_libsrc})
    set_target_properties(libcreativity-gui PROPERTIES
        OUTPUT_NAME creativity-gui
        VERSION "${LIBCRG_CUR}.${LIBCRG_REV}.${LIBCRG_AGE}"
        SOVERSION "${LIBCRG_CUR}"
    )
    foreach(pkg ${creativity_gui_pkgconfig})
        target_link_libraries(libcreativity-gui ${${pkg}_LIBRARIES})
    endforeach()
    list(APPEND creativity_install_targets libcreativity-gui)

    target_link_libraries(libcreativity-gui libcreativity)

    add_executable(creativity-gui gui.cpp)
    target_link_libraries(creativity-gui libcreativity-gui)
    list(APPEND creativity_install_targets creativity-gui)

    # Copy the glade and icon to the install datadir
    install(FILES gui.glade creativity-gui.svg DESTINATION ${CMAKE_INSTALL_DATADIR}/creativity/)

    # Also copy them into the build dir (so that running ./creativity-gui from build/ works)
    configure_file(gui.glade ${CMAKE_BINARY_DIR}/gui.glade COPYONLY)
    configure_file(creativity-gui.svg ${CMAKE_BINARY_DIR}/creativity-gui.svg COPYONLY)

    # Install the desktop file
    install(FILES creativity-gui.desktop DESTINATION ${CMAKE_INSTALL_DATADIR}/applications/)

    # Mime type registration
    install(FILES creativity-gui.xml DESTINATION ${CMAKE_INSTALL_DATADIR}/mime/packages/)

    # Also copy the icon into the pixmaps directory so it'll be picked up by things reading the desktop file
    install(FILES creativity-gui.svg DESTINATION ${CMAKE_INSTALL_DATADIR}/pixmaps/)
    # And put a copy in for the mime type icon
    install(FILES creativity-gui.svg DESTINATION ${CMAKE_INSTALL_DATADIR}/pixmaps/ RENAME application-x-creativity.crstate.svg)

endif()

if (IWYU)
    foreach(tgt ${creativity_install_targets})
        set_property(TARGET ${tgt} PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU_ARGS})
    endforeach()
endif()

if (API_DOCS)
    add_subdirectory(doc)
endif()

add_subdirectory(scratch)

install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/creativity/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/creativity/"
    COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/creativity/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/creativity/"
    COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")

install(TARGETS ${creativity_install_targets}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

install(FILES README.md LICENSE DESTINATION "${CMAKE_INSTALL_DOCDIR}" COMPONENT docs)
