cmake_minimum_required(VERSION 3.11)
project(vash
	VERSION 0.1
	LANGUAGES CXX
)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_EXTENSIONS OFF)
if (NOT DEFINED CMAKE_CXX_STANDARD)
	set(CMAKE_CXX_STANDARD 14)
endif()

# Adding a Profile build type
get_property(isMultiConfig GLOBAL
	PROPERTY GENERATOR_IS_MULTICONFIG
)
if(isMultiConfig)
	if(NOT "Profile" IN_LIST CMAKE_CONFIGURATION_TYPES)
		list(APPEND CMAKE_CONFIGURATION_TYPES Profile)
	endif()
	if(NOT "Test" IN_LIST CMAKE_CONFIGURATION_TYPES)
		list(APPEND CMAKE_CONFIGURATION_TYPES Test)
	endif()
else()
	set(allowedBuildTypes Debug Release Profile Test)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
		STRINGS "${allowedBuildTypes}"
	)
	if(NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE) # default build type is Release
	elseif(NOT CMAKE_BUILD_TYPE IN_LIST allowedBuildTypes)
		message(FATAL_ERROR
			"Unknown build type ${CMAKE_BUILD_TYPE}"
		)
	endif()
endif()
#
# Set warnings
#
set(CLANG_WARNINGS
	-Wall
	-Wextra
	-Wconversion
	-Wpedantic
	-Wshadow
	-Wnon-virtual-dtor
	-Wold-style-cast
	-Wcast-align
	-Woverloaded-virtual
	-Wsign-conversion
	-Wnull-dereference
	-Wdouble-promotion
	-Wformat=2
	-Wimplicit-fallthrough
)
set(GCC_WARNINGS
	${CLANG_WARNINGS} 
	-Wmisleading-indentation
	-Wduplicated-cond
	-Wduplicated-branches
	-Wlogical-op
	-Wuseless-cast
)

if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
	set(PROJECT_WARNINGS_CXX ${CLANG_WARNINGS})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(PROJECT_WARNINGS_CXX ${GCC_WARNINGS})
else()
	message(AUTHOR_WARNING "No compiler warnings set for CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()
#
# Tests with Catch2
#
option(BUILD_TESTS
	"Build tests"
	OFF
)
if(BUILD_TESTS)
	Include(FetchContent)

	FetchContent_Declare(
	  Catch2
	  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	  GIT_TAG        v3.3.2 # or a later release
	)
	FetchContent_MakeAvailable(Catch2)
	list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()
#
# Find active available sanitizers
#
set(SANITIZER_LIST "")
if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(SANITIZERS "")
	check_cxx_compiler_flag(-fsanitize=address   HAS_ADDRESS_SANITIZER)
	check_cxx_compiler_flag(-fsanitize=leak      HAS_LEAK_SANITIZER)
	check_cxx_compiler_flag(-fsanitize=undefined HAS_UNDEFINED_SANITIZER)
	if(${HAS_ADDRESS_SANITIZER})
		list(APPEND SANITIZERS "address")
	endif()
	if(${HAS_LEAK_SANITIZER})
		list(APPEND SANITIZERS "leak")
	endif()
	if(${HAS_UNDEFINED_SANITIZER})
		list(APPEND SANITIZERS "undefined")
	endif()
	list(JOIN
		SANITIZERS
		","
		SANITIZER_LIST)
endif()


set(CMAKE_CXX_FLAGS_PROFILE
	"-g -p -O3"
)
set(CMAKE_CXX_SHARED_LINKER_FLAGS_PROFILE
	"-g -p -O3"
)
set(CMAKE_CXX_STATIC_LINKER_FLAGS_PROFILE
	"-g -p -O3"
)

set(CMAKE_CXX_FLAGS_DEBUG
	"-g -O1 -fsanitize=${SANITIZER_LIST}"
)
set(CMAKE_CXX_SHARED_LINKER_FLAGS_DEBUG
	"-g -O1 -fsanitize=${SANITIZER_LIST}"
)
set(CMAKE_CXX_STATIC_LINKER_FLAGS_DEBUG
	"-g -O1 -fsanitize=${SANITIZER_LIST}"
)

set(CMAKE_CXX_FLAGS_TEST
	"-g -O3 -fsanitize=${SANITIZER_LIST}"
)
set(CMAKE_CXX_SHARED_LINKER_FLAGS_TEST
	"-g -O3 -fsanitize=${SANITIZER_LIST}"
)
set(CMAKE_CXX_STATIC_LINKER_FLAGS_TEST
	"-g -O3 -fsanitize=${SANITIZER_LIST}"
)

include(GNUInstallDirs)
# executable
add_executable(ldblocks
	apps/ldblocks.cpp
)
target_include_directories(ldblocks
	PRIVATE externals/bayesicUtilities/include
	PRIVATE include
)
target_link_libraries(ldblocks
	PRIVATE vash
)
target_compile_options(ldblocks
	PRIVATE ${PROJECT_WARNINGS_CXX}
)
set_target_properties(ldblocks PROPERTIES
	CXX_STANDARD_REQUIRED ON
)
install(TARGETS ldblocks
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# library
add_library(vash
	src/gvarHash.cpp
	src/vashFunctions.cpp
	externals/bayesicUtilities/src/random.cpp
)

target_include_directories(vash
	PRIVATE externals/bayesicUtilities/include
	PRIVATE include
)
set_target_properties(vash PROPERTIES
	POSITION_INDEPENDENT_CODE ON
)
target_compile_options(vash
	PRIVATE ${PROJECT_WARNINGS_CXX} -mbmi -mbmi2 -mpopcnt -pthread
)
set_target_properties(vash PROPERTIES
	CXX_STANDARD_REQUIRED ON
)
# no install for now

# tests
if(BUILD_TESTS)
	add_executable(tests
		tests/tests.cpp
	)
	target_link_libraries(tests
		PRIVATE vash
		PRIVATE Catch2::Catch2WithMain
	)
	target_include_directories(tests
		PRIVATE include
		PRIVATE externals/bayesicUtilities/include
	)
	include(CTest)
	include(Catch)
	catch_discover_tests(tests)
endif()
