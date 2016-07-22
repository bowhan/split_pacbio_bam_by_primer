cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0048 NEW)
project(SPLIT_PRIMER_FROM_PBBAM VERSION 0.0.0.1 LANGUAGES CXX)

string(TOLOWER ${CMAKE_PROJECT_NAME} MAIN_EXE_NAME)

# variables to all sorts of directories
set(${CMAKE_PROJECT_NAME}_ROOT_DIR ${PROJECT_SOURCE_DIR})
set(INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)
set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/third-party)

# configure files
configure_file(${INCLUDE_DIRS}/config.h.in ${INCLUDE_DIRS}/config.h)

# CXX flags
set(${CMAKE_PROJECT_NAME}_FLAGS "-std=c++11 -Wall -Wextra -Werror=strict-aliasing -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${${CMAKE_PROJECT_NAME}_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${${CMAKE_PROJECT_NAME}_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS_RELEASE} ${${CMAKE_PROJECT_NAME}_FLAGS}")
if (APPLE)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
else ()
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
endif ()

#  pbcopper
add_subdirectory(${THIRD_PARTY_DIR}/pbcopper external/pbcopper/build)

find_package(htslib
        PATHS ${THIRD_PARTY_DIR}/htslib
        REQUIRED)

# pbbam
# prevent pbbam from building unnecessaries
set(PacBioBAM_build_docs OFF CACHE INTERNAL "" FORCE)
set(PacBioBAM_build_tests OFF CACHE INTERNAL "" FORCE)
set(PacBioBAM_build_tools OFF CACHE INTERNAL "" FORCE)
add_subdirectory(${THIRD_PARTY_DIR}/pbbam external/pbbam/build)

# complete striped SWmake in
add_subdirectory("${THIRD_PARTY_DIR}/Complete-Striped-Smith-Waterman-Library" external/ssw/build)

# BOOST
find_package(Boost REQUIRED)

# include
include_directories(
        ${INCLUDE_DIRS}
        ${THIRD_PARTY_DIR}
        ${Boost_INCLUDE_DIRS}
)

# exe
add_executable(${MAIN_EXE_NAME}
        ${SOURCE_DIR}/main.cpp
        ${SOURCE_DIR}/common.cpp
        )

# link
target_link_libraries(${MAIN_EXE_NAME}
        pbcopper
        pbbam
        sswlib
        ${ZLIB_LIBRARIES}
        )
add_dependencies(${MAIN_EXE_NAME} htslib)
set_target_properties(${MAIN_EXE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
