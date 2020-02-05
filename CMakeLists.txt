cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# version
foreach (level MAJOR MINOR PATCH)
    file(STRINGS include/version.inc
            _define_${level}
            REGEX "#define PROG_${level}_VERSION")
    string(REGEX MATCH "([0-9]+)" _${level} "${_define_${level}}")
endforeach ()

set(VERSION "${_MAJOR}.${_MINOR}.${_PATCH}")
project(SPLIT_PRIMER_FROM_PBBAM VERSION ${VERSION} LANGUAGES CXX)
string(TOLOWER ${CMAKE_PROJECT_NAME} MAIN_EXE_NAME)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(${CMAKE_PROJECT_NAME}_ROOT_DIR ${PROJECT_SOURCE_DIR})
set(INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)
set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/third-party)

find_package(Threads REQUIRED)
find_package(Intellia REQUIRED)

SET(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS filesystem REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED IMPORTED_TARGET "htslib")
pkg_check_modules(PBBAM REQUIRED IMPORTED_TARGET "pbbam")


add_subdirectory("${THIRD_PARTY_DIR}/Complete-Striped-Smith-Waterman-Library" external/ssw/build)

# exe
add_executable(${MAIN_EXE_NAME}
        ${SOURCE_DIR}/main.cpp
        ${SOURCE_DIR}/common.cpp
        )
#
target_include_directories(${MAIN_EXE_NAME}
        PUBLIC
        ${INCLUDE_DIRS}
        ${HTSLIB_INCLUDEDIR}
        ${PBBAM_INCLUDEDIR}
        ${THIRD_PARTY_DIR}
        )


# link
target_link_directories(${MAIN_EXE_NAME}
        PUBLIC
        ${HTSLIB_LIBRARY_DIRS}
        ${PBBAM_LIBRARY_DIRS}
        )

target_link_libraries(${MAIN_EXE_NAME}
        ${HTSLIB_LIBRARIES}
        ${PBBAM_LIBRARIES}
        ${Boost_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        sswlib
        )

# install
install(TARGETS ${MAIN_EXE_NAME}
        RUNTIME DESTINATION bin
        )