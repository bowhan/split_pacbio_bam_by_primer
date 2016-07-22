cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0048 NEW)
project(SPLIT_POLYT_FROM_PBBAM VERSION 0.0.0.1 LANGUAGES CXX)

string(TOLOWER ${CMAKE_PROJECT_NAME} MAIN_EXE_NAME)

# variables to all sorts of directories
set(${CMAKE_PROJECT_NAME}_ROOT_DIR ${PROJECT_SOURCE_DIR})
set(${CMAKE_PROJECT_NAME}_INCLUDE_DIRS ${${CMAKE_PROJECT_NAME}_ROOT_DIR}/include)
set(${CMAKE_PROJECT_NAME}_SOURCE_DIR ${${CMAKE_PROJECT_NAME}_ROOT_DIR}/src)
set(${CMAKE_PROJECT_NAME}_THIRD_PARTY_DIR ${${CMAKE_PROJECT_NAME}_ROOT_DIR}/third-party)
set(${CMAKE_PROJECT_NAME}_TEST_DIR ${${CMAKE_PROJECT_NAME}_ROOT_DIR}/test)

# configure files
configure_file(${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}/config.h.in ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}/config.h)

# CXX flags
set(${CMAKE_PROJECT_NAME}_FLAGS " -static -std=c++11 -Wall -Wextra -Werror=strict-aliasing -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${${CMAKE_PROJECT_NAME}_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${${CMAKE_PROJECT_NAME}_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-static ${CMAKE_CXX_FLAGS_RELEASE} ${${CMAKE_PROJECT_NAME}_FLAGS}")
if (APPLE)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
else ()
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
endif ()

#  COMMON INTERFACE
add_library(COMMON_INTERFACE INTERFACE)

#  pbcopper
add_subdirectory(${${CMAKE_PROJECT_NAME}_THIRD_PARTY_DIR}/pbcopper external/pbcopper/build)
target_link_libraries(COMMON_INTERFACE INTERFACE pbcopper)

# pbbam
# prevent pbbam from building unnecessaries
set(PacBioBAM_build_docs OFF CACHE INTERNAL "" FORCE)
set(PacBioBAM_build_tests OFF CACHE INTERNAL "" FORCE)
set(PacBioBAM_build_tools OFF CACHE INTERNAL "" FORCE)
add_subdirectory(${${CMAKE_PROJECT_NAME}_THIRD_PARTY_DIR}/pbbam external/pbbam/build)
target_link_libraries(COMMON_INTERFACE INTERFACE pbbam)

# complete striped SW
add_subdirectory("third-party/Complete-Striped-Smith-Waterman-Library" external/ssw/build)
target_link_libraries(COMMON_INTERFACE INTERFACE sswlib)
# BOOST
find_package(Boost REQUIRED)

# include
include_directories(
        ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
        ${${CMAKE_PROJECT_NAME}_THIRD_PARTY_DIR}
        ${Boost_INCLUDE_DIRS}
)

# exe
add_executable(${MAIN_EXE_NAME}
        ${${CMAKE_PROJECT_NAME}_SOURCE_DIR}/main.cpp
        ${${CMAKE_PROJECT_NAME}_SOURCE_DIR}/common.cpp
        )

# link
target_link_libraries(${MAIN_EXE_NAME}
        PUBLIC
        COMMON_INTERFACE
        ${ZLIB_LIBRARIES}
        )

set(ALL_EXE_TARGETS ${MAIN_EXE_NAME})

install(TARGETS ${ALL_EXE_TARGETS}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)