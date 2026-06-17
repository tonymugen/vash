find_package(Doxygen REQUIRED)

file(GLOB PUBLIC_HEADERS ${PROJECT_SOURCE_DIR}/plugin/*.hpp)

set(DOXYGEN_INPUT_DIR  ${PROJECT_SOURCE_DIR}/plugin)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doxygen)
set(README_FILE        ${PROJECT_SOURCE_DIR}/README.md)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN        ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT       ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
# Custom stylesheet matching the bayesicresearch.org design (Nord palette + the
# site's typefaces). Self-contained, so the build needs no external theme.
set(CUSTOM_CSS_FILE    ${CMAKE_CURRENT_SOURCE_DIR}/bayesic.css)

# Replace variables in Doxyfile
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
add_custom_command(
	OUTPUT            ${DOXYGEN_INDEX_FILE}
	DEPENDS           ${PUBLIC_HEADERS}
	                  ${README_FILE}
	                  ${CUSTOM_CSS_FILE}
	COMMAND           ${DOXYGEN_EXECUTABLE}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	MAIN_DEPENDENCY   ${DOXYFILE_IN}
	COMMENT           "Generating docs"
	VERBATIM
)
add_custom_target(Doxygen
	ALL
	DEPENDS ${DOXYGEN_INDEX_FILE}
)
