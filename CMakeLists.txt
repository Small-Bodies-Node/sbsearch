cmake_minimum_required(VERSION 3.10)
project(sbsearch
    VERSION 3.0.0)
configure_file(src/sbsearch.h.in sbsearch.h)
include(GoogleTest)
include(CTest)

# compile to ephemeris.o rather than ephemeris.cc.o
# this will help gcov find the compiled files
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)

add_executable(build_test_db
    src/ephemeris.cc
    src/build_test_db.cc
    src/observation.cc
    src/sbsearch_testing.cc
    src/sbsdb.cc
    src/sbsdb_sqlite3.cc
    src/util.cc
    )
target_include_directories(build_test_db PUBLIC
"${PROJECT_BINARY_DIR}"
)
target_link_libraries(build_test_db s2 m sqlite3)

add_executable(query_test_db
    src/ephemeris.cc
    src/observation.cc
    src/query_test_db.cc
    src/sbsearch_testing.cc
    src/sbsdb.cc
    src/sbsdb_sqlite3.cc
    src/util.cc
)
target_include_directories(query_test_db PUBLIC
"${PROJECT_BINARY_DIR}"
)
target_link_libraries(query_test_db s2 m sqlite3 s2testing)

add_executable(sbsearch_test
    src/ephemeris.cc
    src/ephemeris_test.cc
    src/observation.cc
    src/observation_test.cc
    src/sbsearch_testing.cc
    src/test.cc
    src/util.cc
    src/util_test.cc
)
target_include_directories(sbsearch_test PUBLIC
"${PROJECT_BINARY_DIR}"
)
target_compile_options(sbsearch_test PUBLIC
    -g
    -O0
    # -Wall
    --coverage
)
target_link_libraries(sbsearch_test s2 m sqlite3 gtest -lgcov)

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
