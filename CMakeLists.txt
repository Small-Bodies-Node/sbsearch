cmake_minimum_required(VERSION 3.10)
project(sbsearch
    VERSION 4.0.0)
configure_file(${CMAKE_SOURCE_DIR}/src/config.h.in ${CMAKE_BINARY_DIR}/config.h)
include(GoogleTest)
include(CTest)
find_package(Boost 1.81.0 REQUIRED COMPONENTS program_options json filesystem)

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

find_library(S2LIB NAME s2
    HINTS "${CMAKE_PREFIX_PATH}/lib"
)
if (NOT S2LIB)
    message(FATAL_ERROR "Cannot find libs2")
endif()
message(STATUS "Found libs2 (${S2LIB})")

find_path(S2INCLUDE NAME s2/s2polygon.h
    HINTS "${CMAKE_PREFIX_PATH}"
)
message(STATUS "Found s2 include directory (${S2INCLUDE})")

include_directories(
    src/
    ${S2INCLUDE}
    ${PROJECT_BINARY_DIR}
    ${CMAKE_BINARY_DIR}
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# handy to avoid excessive warnings due to s2 version
add_definitions(-Wno-deprecated-declarations)

# compile to ephemeris.o rather than ephemeris.cc.o
# this will help gcov find the compiled files
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)

add_library(sbsearch SHARED
    src/cli.cc
    src/csv.cc
    src/date.cc
    src/ephemeris.cc
    src/env.cc
    src/files.cc
    src/found.cc
    src/horizons.cc
    src/intersection.cc
    src/indexer.cc
    src/logging.cc
    src/moving_target.cc
    src/observation.cc
    src/query_info.cc
    src/sbsdb/add.cc
    src/sbsdb/count.cc
    src/sbsdb/find.cc
    src/sbsdb/get.cc
    src/sbsdb/postgresql.cc
    src/sbsdb/remove.cc
    src/sbsdb/update.cc
    src/sbsdb/verify.cc
    src/sbsearch/add.cc
    src/sbsearch/find_fixed_target.cc
    src/sbsearch/find_moving_target.cc
    src/sbsearch/index.cc
    src/sbsearch/sbsearch.cc
    src/sbsearch/update.cc
    src/table.cc
    src/util/math.cc
    src/util/polygon.cc
    src/util/spherical.cc
    src/util/string.cc
    src/sofa/anp.c
    src/sofa/cal2jd.c
    src/sofa/dat.c
    src/sofa/era00.c
    src/sofa/d2tf.c
    src/sofa/jd2cal.c
    src/sofa/taiut1.c
    src/sofa/taiutc.c
    src/sofa/utctai.c
    src/sofa/utcut1.c
)

set(SBSEARCHLIBS
    sbsearch
    Boost::program_options
    Boost::filesystem
    Boost::json
    curl
    m
    ${S2LIB}
    sqlite3
    pq
    pqxx
)

add_executable(sbs-ephemeris src/sbs-ephemeris.cc)
target_link_libraries(sbs-ephemeris
    ${SBSEARCHLIBS}
)

add_executable(sbs-configure src/sbs-configure.cc)
target_link_libraries(sbs-configure
    ${SBSEARCHLIBS}
)

add_executable(sbs-found src/sbs-found.cc)
target_link_libraries(sbs-found
    ${SBSEARCHLIBS}
)

add_executable(sbs-moving-target src/sbs-moving-target.cc)
target_link_libraries(sbs-moving-target
    ${SBSEARCHLIBS}
)

add_executable(sbs-observation
    src/sbs-observation/add.cc
    src/sbs-observation/arguments.cc
    src/sbs-observation/main.cc
    src/sbs-observation/summary.cc
)
target_link_libraries(sbs-observation
    ${SBSEARCHLIBS}
)

add_executable(sbs-observatory src/sbs-observatory.cc)
target_link_libraries(sbs-observatory
    ${SBSEARCHLIBS}
)
    
add_executable(sbs-query src/sbs-query.cc)
target_link_libraries(sbs-query
    ${SBSEARCHLIBS}
)
    
add_executable(build_test_db src/tests/build_test_db.cc)
target_link_libraries(build_test_db PUBLIC
    ${SBSEARCHLIBS}
)

add_executable(query_test_db src/tests/query_test_db.cc)
target_link_libraries(query_test_db PUBLIC
    ${SBSEARCHLIBS}
)

add_executable(sbsearch_test
    src/tests/csv_test.cc
    src/tests/date_test.cc
    src/tests/ephemeris_test.cc
    src/tests/files_test.cc
    src/tests/found_test.cc
    src/tests/horizons_test.cc
    src/tests/intersection_test.cc
    src/tests/indexer_test.cc
    src/tests/logging_test.cc
    src/tests/moving_target_test.cc
    src/tests/observation_test.cc
    src/tests/observatory_test.cc
    src/tests/sbsearch_test.cc
    src/tests/sbsdb/ephemeris_io_test.cc
    src/tests/sbsdb/find_test.cc
    src/tests/sbsdb/found_io_test.cc
    src/tests/sbsdb/indexer_options_io_test.cc
    src/tests/sbsdb/moving_target_io_test.cc
    src/tests/sbsdb/observations_io_test.cc
    src/tests/sbsdb/observatory_io_test.cc
    src/tests/sbsdb/postgresql_test.cc
    src/tests/table_test.cc
    src/tests/test.cc
    src/tests/util/math_test.cc
    src/tests/util/polygon_test.cc
    src/tests/util/spherical_test.cc
    src/tests/util/string_test.cc
)
target_compile_options(sbsearch_test PUBLIC
    -g
    -O0
    # -Wall
    --coverage
)
target_link_libraries(sbsearch_test
    ${SBSEARCHLIBS}
    gtest
    -lgcov
)

add_executable(diagnostics
    src/diagnostics.cc
)
target_link_libraries(diagnostics
    ${SBSEARCHLIBS}
)

install(TARGETS
    sbsearch
    sbs-ephemeris
    sbs-configure
    sbs-found
    sbs-moving-target
    sbs-observation
    sbs-observatory
    sbs-query
)

############# Examples
add_custom_target(examples
    DEPENDS sbs-add-lsst
)

add_executable(sbs-add-lsst src/sbs-add-lsst.cc)
target_link_libraries(sbs-add-lsst
    ${SBSEARCHLIBS}
)

############# Coverage targets
# Based on https://dr-kino.github.io/2019/12/22/test-coverage-using-gtest-gcov-and-lcov/

set(OBJECT_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/sbsearch_test.dir/src)

add_custom_target(gcov
    COMMAND mkdir -p gcoverage
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS sbsearch_test
)
add_custom_command(TARGET gcov
    COMMAND sbsearch_test
    COMMAND echo "=================== GCOV ===================="
    COMMAND gcov -b ${CMAKE_SOURCE_DIR}/src/*.cc -o ${OBJECT_DIR}
    COMMAND echo "-- Source diretories: ${CMAKE_SOURCE_DIR}/src/"
    COMMAND echo "-- Coverage files have been output to ${CMAKE_BINARY_DIR}/gcoverage"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gcoverage
)
add_dependencies(gcov sbsearch_test)

add_custom_target(lcov
    COMMAND mkdir -p lcoverage
    DEPENDS gcov
)
add_custom_command(TARGET lcov
    COMMAND echo "=================== LCOV ===================="
    COMMAND echo "-- Passing lcov tool under code coverage"
    COMMAND lcov --capture --directory ../ --output-file lcoverage/main_coverage.info
    COMMAND echo "-- Generating HTML output files"
    COMMAND genhtml lcoverage/main_coverage.info --output-directory lcoverage
)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gcoverage)

add_custom_target(cov-clean
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    COMMAND rm -f ${OBJECT_DIR}/*.gcno
    COMMAND rm -f ${OBJECT_DIR}/*.gcda
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Testing
enable_testing()
gtest_discover_tests(sbsearch_test)
