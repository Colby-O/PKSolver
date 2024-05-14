message(STATUS "PARAMETER_EXTRACTION - Lading Common Configuration...")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

include_directories("${ROOT_DIR}/include/")

file(GLOB_RECURSE COMMON_DIR_HEADERS "${CMAKE_CURRENT_LIST_DIR}/*.hpp")
file(GLOB_RECURSE COMMON_DIR_HEADERS "${CMAKE_CURRENT_LIST_DIR}/*.h")
add_custom_target(common_headers SOURCES ${COMMON_DIR_HEADERS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lstdc++ -std=c++17 -O0 -g")
