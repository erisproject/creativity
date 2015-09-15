
file(GLOB scratch_progs_src RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.cpp)

string(REGEX REPLACE ".cpp(;|$)" "\\1" scratch_progs "${scratch_progs_src}")

include_directories("${CMAKE_SOURCE_DIR}")

foreach(prog ${scratch_progs})
    add_executable(${prog} EXCLUDE_FROM_ALL "${prog}.cpp")

    foreach(pkg ERIS)
        target_link_libraries(${prog} ${${pkg}_LIBRARIES})
    endforeach()
    target_link_libraries(${prog} libcreativity)
    target_link_libraries(${prog} libcreativity-data)
endforeach()

add_custom_target(scratch DEPENDS ${scratch_progs})

