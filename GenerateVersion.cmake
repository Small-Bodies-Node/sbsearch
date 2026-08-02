# h/t https://www.marcusfolkesson.se/blog/git-version-in-cmake/
find_package(Git)

if(GIT_EXECUTABLE)
  get_filename_component(WORKING_DIR ${SRC} DIRECTORY)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
    WORKING_DIRECTORY ${WORKING_DIR}
    OUTPUT_VARIABLE SBSEARCH_VERSION
    RESULT_VARIABLE ERROR_CODE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()

if(SBSEARCH_VERSION STREQUAL "")
  set(SBSEARCH_VERSION 0.0.0-unknown)
  message(WARNING "Failed to determine version from Git tags. Using default version \"${SBSEARCH_VERSION}\".")
endif()

configure_file(${SRC} ${DST} @ONLY)
