cmake_minimum_required(VERSION 3.10)
project(sbsearch
    VERSION 3.0.0)
configure_file(src/sbsearch.h.in sbsearch.h)

add_executable(build_test_db
    src/sbsearch.h
    src/sbsdb_sqlite3.cc
    src/util.h
    src/util.cc
    src/build_test_db.cc
    src/observation.h
    src/observation.cc
    src/sbsdb.h
    src/sbsdb.cc
    src/sbsdb_sqlite3.h
)
target_include_directories(build_test_db PUBLIC
"${PROJECT_BINARY_DIR}"
)
target_link_libraries(build_test_db s2 m sqlite3)

add_executable(query_test_db
    src/query_test_db.cc
    src/sbsearch.h
    src/util.h
    src/util.cc
    src/observation.h
    src/observation.cc
    src/ephemeris.h
    src/ephemeris.cc
    src/sbsearch_testing.h
    src/sbsearch_testing.cc
)
target_include_directories(query_test_db PUBLIC
"${PROJECT_BINARY_DIR}"
)
target_link_libraries(query_test_db s2 m sqlite3 s2testing)

add_executable(test
    src/test.cc
    src/sbsearch.h
    src/util.h
    src/util.cc
    src/observation.h
    src/observation.cc
    src/ephemeris.h
    src/ephemeris.cc
    src/sbsearch_testing.h
    src/sbsearch_testing.cc
    src/ephemeris_test.cc
    src/observation_test.cc
)
target_include_directories(test PUBLIC
"${PROJECT_BINARY_DIR}"
)
target_link_libraries(test s2 m sqlite3 gtest)

# install(TARGETS sbsearch DESTINATION bin)
# install(FILES "${PROJECT_BINARY_DIR}/sbsearch.h"
#   DESTINATION include
#   )
