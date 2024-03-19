cmake_minimum_required(VERSION 3.10)
project(sbsearch
    VERSION 3.0.0)
configure_file(src/config.h.in config.h)
include(GoogleTest)
include(CTest)
find_package(Boost 1.81.0 REQUIRED COMPONENTS program_options json)

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
    "${S2INCLUDE}"
    "${PROJECT_BINARY_DIR}"
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# compile to ephemeris.o rather than ephemeris.cc.o
# this will help gcov find the compiled files
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)

add_library(sbsearch SHARED
    src/ephemeris.cc
    src/indexer.cc
    src/logging.cc
    src/moving_target.cc
    src/observation.cc
    src/sbs-cli.cc
    src/sbsdb.cc
    src/sbsdb_sqlite3.cc
    src/sbsearch.cc
    src/util.cc
    src/sofa/anp.c
    src/sofa/cal2jd.c
    src/sofa/dat.c
    src/sofa/era00.c
    src/sofa/jd2cal.c
    src/sofa/taiut1.c
    src/sofa/taiutc.c
    src/sofa/utctai.c
    src/sofa/utcut1.c
)

add_executable(sbs-ephemeris src/sbs-ephemeris.cc)
target_link_libraries(sbs-ephemeris
    sbsearch
    boost_program_options
    curl
    m
    ${S2LIB}
    sqlite3
)

add_executable(sbs-configure src/sbs-configure.cc)
target_link_libraries(sbs-configure
    sbsearch
    boost_program_options
    curl
    m
    ${S2LIB}
    sqlite3
)

add_executable(sbs-moving-target src/sbs-moving-target.cc)
target_link_libraries(sbs-moving-target
    sbsearch
    boost_program_options
    curl
    m
    ${S2LIB}
    sqlite3
)

add_executable(sbs-observation src/sbs-observation.cc)
target_link_libraries(sbs-observation
    sbsearch
    boost_program_options
    curl
    m
    ${S2LIB}
    sqlite3
)

add_executable(sbs-observatory src/sbs-observatory.cc)
target_link_libraries(sbs-observatory
    sbsearch
    boost_program_options
    curl
    m
    ${S2LIB}
    sqlite3
)
    
add_executable(sbs-query src/sbs-query.cc)
target_link_libraries(sbs-query
    sbsearch
    boost_program_options
    curl
    m
    ${S2LIB}
    sqlite3
)
    
add_executable(build_test_db src/build_test_db.cc)
target_link_libraries(build_test_db PUBLIC
    sbsearch
    boost_program_options
    curl
    m
    ${S2LIB}
    sqlite3
)

add_executable(query_test_db src/query_test_db.cc)
target_link_libraries(query_test_db PUBLIC
    sbsearch
    boost_program_options
    curl
    m
    ${S2LIB}
    sqlite3
)

add_executable(reindex_test_db src/reindex_test_db.cc)
target_link_libraries(reindex_test_db
    sbsearch
    boost_program_options
    curl
    m
    ${S2LIB}
    sqlite3
)

add_executable(sbsearch_test
    src/ephemeris_test.cc
    src/indexer_test.cc
    src/logging_test.cc
    src/moving_target_test.cc
    src/observation_test.cc
    src/observatory_test.cc
    src/sbsearch_test.cc
    src/sbsdb_sqlite3_test.cc
    src/test.cc
    src/util_test.cc
)
target_compile_options(sbsearch_test PUBLIC
    -g
    -O0
    # -Wall
    --coverage
)
target_link_libraries(sbsearch_test
    sbsearch
    boost_program_options
    curl
    ${S2LIB}
    m
    sqlite3
    gtest
    -lgcov
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
