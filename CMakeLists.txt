cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# version
foreach (level MAJOR MINOR PATCH)
    file(STRINGS include/version.inc
            _define_${level}
            REGEX "#define PROG_${level}_VERSION")
    string(REGEX MATCH "([0-9]+)" _${level} "${_define_${level}}")
endforeach ()

set(VERSION "${_MAJOR}.${_MINOR}.${_PATCH}")

# project
project(SPLIT_PRIMER_FROM_PBBAM VERSION ${VERSION} LANGUAGES C CXX)
string(TOLOWER ${CMAKE_PROJECT_NAME} MAIN_EXE_NAME)

# flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# project-specific variables
set(${CMAKE_PROJECT_NAME}_ROOT_DIR ${PROJECT_SOURCE_DIR})
set(INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)
set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/third-party)

# required packages
# 1. threads
find_package(Threads REQUIRED)
# 2. boost
SET(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS filesystem REQUIRED)
# 3. search for htslib and pbbam with pkg-config
find_package(PkgConfig REQUIRED)
pkg_search_module(HTS REQUIRED IMPORTED_TARGET "htslib")
pkg_search_module(PBBAM REQUIRED IMPORTED_TARGET "pbbam")

# exe
add_executable(${MAIN_EXE_NAME}
        ${SOURCE_DIR}/main.cpp
        ${SOURCE_DIR}/common.cpp
        ${SOURCE_DIR}/impl/ssw/ssw_impl.c
        ${SOURCE_DIR}/Ssw.cpp
        )
# include
target_include_directories(${MAIN_EXE_NAME}
        PUBLIC
        ${INCLUDE_DIRS}
        ${HTSLIB_INCLUDEDIR}
        ${PBBAM_INCLUDEDIR}
        ${THIRD_PARTY_DIR}
        )
# linking
target_link_libraries(${MAIN_EXE_NAME}
        PkgConfig::HTS
        PkgConfig::PBBAM
        ${Boost_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        )
# install
install(TARGETS ${MAIN_EXE_NAME}
        RUNTIME DESTINATION bin
        )