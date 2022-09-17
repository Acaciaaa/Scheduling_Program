cmake_minimum_required (VERSION 3.5)

project (hls)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
include_directories(${PROJECT_SOURCE_DIR}/src)

include_directories(${PROJECT_SOURCE_DIR}/glpk/include)
find_library(GLPK_LIB glpk HINTS ${PROJECT_SOURCE_DIR}/glpk/lib)

add_executable(hls ${SRC_LIST})

target_link_libraries(hls ${GLPK_LIB})