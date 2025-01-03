cmake_minimum_required(VERSION 3.16)
project(BlackScholesPDE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(${CMAKE_SOURCE_DIR}/src)

file(GLOB SRC_FILES
    ${CMAKE_SOURCE_DIR}/src/*.cpp
)

file(GLOB HEADER_FILES
    ${CMAKE_SOURCE_DIR}/src/*.hpp
)

add_executable(BlackScholesPDE ${SRC_FILES})
set_target_properties(BlackScholesPDE PROPERTIES OUTPUT_NAME main)

file(GLOB EXAMPLE_FILES
    ${CMAKE_SOURCE_DIR}/examples/*.cpp
)

foreach(EXAMPLE_FILE ${EXAMPLE_FILES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_FILE} ${SRC_FILES})
endforeach()

file(GLOB TEST_FILES
    ${CMAKE_SOURCE_DIR}/tests/*.cpp
)

foreach(TEST_FILE ${TEST_FILES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE} ${SRC_FILES})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

enable_testing()
