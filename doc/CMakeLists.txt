
# Need at least version 1.8.2 (as that was when c++11 support was added to doxygen)
find_package(Doxygen 1.8.2 REQUIRED)

if(NOT DOXYGEN_DOT_FOUND)
    message(FATAL_ERROR "Could NOT find dot (required for doxygen documentation generation)")
endif()

# Can set MATHJAX_CDN to use the MathJax CDN (instead of requiring a local
# MathJax).
set(mathjax_usable_path "")
if (MATHJAX_CDN)
    message(STATUS "Using public MathJax CDN")
else ()
    find_package(MathJax 2.3 REQUIRED)
    set(mathjax_usable_path "MATHJAX_RELPATH = ${MATHJAX_PATH}")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml.in ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml @ONLY)

# Use a custom command to output doxygen's generated index file.  Doing it this way ensures
# that 'make api_doc' only runs doxygen when something in src or include has changed.
add_custom_command(OUTPUT "api/index.html"
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    DEPENDS "Doxyfile.in" "DoxygenLayout.xml.in" ${creativity_doc_include}
        ${CMAKE_SOURCE_DIR}/README.md
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Generating API documentation (using doxygen)"
    )

add_custom_target(api_doc ALL DEPENDS "api/index.html")

set_property(DIRECTORY APPEND PROPERTY
    ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_CURRENT_BINARY_DIR}/api")



install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/api"
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
    COMPONENT docs)
