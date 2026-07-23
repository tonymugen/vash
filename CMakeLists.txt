cmake_minimum_required(VERSION 3.21)
project(vash
	VERSION 0.3
	LANGUAGES CXX
)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_EXTENSIONS OFF)
if (NOT DEFINED CMAKE_CXX_STANDARD)
	set(CMAKE_CXX_STANDARD 17)
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
Include(FetchContent)
set (FETCHCONTENT_QUIET FALSE)
FetchContent_Declare(
	brutilities
	GIT_REPOSITORY https://github.com/tonymugen/bayesicUtilities
	GIT_TAG main
	GIT_PROGRESS TRUE
	SYSTEM            # to disable warnings
)
FetchContent_MakeAvailable(brutilities)

#
# Tests with Catch2
#
option(BUILD_TESTS
	"Build tests"
	OFF
)
if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTS)
	FetchContent_Declare(
		Catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG        v3.5.2 # or a later release
		GIT_PROGRESS   TRUE
		SYSTEM                # to disable Catch2 warnings
	)
	FetchContent_MakeAvailable(Catch2)
endif()

#
# Benchmarks (benchmarking branch only): standalone timing harnesses. Enabling this defines
# VASH_BENCHMARK on the vash library so the in-source [vash-bench] phase markers become active.
#
option(BUILD_BENCHMARKS
	"Build benchmark harnesses and activate in-source phase timers (VASH_BENCHMARK)"
	OFF
)

#
# Documentation
#
option(BUILD_DOCS
	"Build docs"
	OFF
)
if(PROJECT_IS_TOP_LEVEL AND BUILD_DOCS)
	add_subdirectory(docs)
endif()
#
# Find active available sanitizers
#
set(SANITIZER_LIST "")
if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(SANITIZERS "")
	# The sanitizer runtimes are resolved at link time, so the feature probe must link with the
	# flag too: otherwise -fsanitize=address fails the check (unresolved __asan_* references) while
	# compile-only leak spuriously passes. Clear the cached results so the checks re-run with the
	# corrected link options even in an existing build tree.
	unset(HAS_ADDRESS_SANITIZER CACHE)
	unset(HAS_LEAK_SANITIZER CACHE)
	unset(HAS_UNDEFINED_SANITIZER CACHE)
	set(_vash_saved_required_link_options "${CMAKE_REQUIRED_LINK_OPTIONS}")
	set(CMAKE_REQUIRED_LINK_OPTIONS -fsanitize=address)
	check_cxx_compiler_flag(-fsanitize=address    HAS_ADDRESS_SANITIZER)
	set(CMAKE_REQUIRED_LINK_OPTIONS -fsanitize=leak)
	check_cxx_compiler_flag(-fsanitize=leak       HAS_LEAK_SANITIZER)
	# UBSan also resolves its handlers (__ubsan_handle_*) at link time; the earlier "weird linking
	# errors with UBSAN" were the compile-only setup leaving those unresolved, fixed by the probe
	# (and target) now linking with the flag.
	set(CMAKE_REQUIRED_LINK_OPTIONS -fsanitize=undefined)
	check_cxx_compiler_flag(-fsanitize=undefined  HAS_UNDEFINED_SANITIZER)
	set(CMAKE_REQUIRED_LINK_OPTIONS "${_vash_saved_required_link_options}")
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
	# Instrumented translation units compile with these flags. -fno-sanitize-recover=undefined makes
	# UBSan trap (abort) on the first error instead of printing and continuing, so undefined behavior
	# fails the test run rather than slipping past a green suite. It is a codegen flag (compile only);
	# the link step needs just -fsanitize=${SANITIZER_LIST} to pull in the runtimes.
	set(SANITIZER_COMPILE_FLAGS -fsanitize=${SANITIZER_LIST})
	if(HAS_UNDEFINED_SANITIZER)
		list(APPEND SANITIZER_COMPILE_FLAGS -fno-sanitize-recover=undefined)
	endif()
endif()

set(CMAKE_CXX_FLAGS_RELEASE
	"-O3 -DNDEBUG"
)

set(CMAKE_CXX_FLAGS_PROFILE
	"-g -p -O3 -DNDEBUG"
)

set(CMAKE_CXX_FLAGS_DEBUG
	"-g -O1"
)

set(CMAKE_CXX_FLAGS_TEST
	"-g -O3"
)

#
# Parallel-algorithm (PSTL) backend
#
# std::execution::par needs a backend, and only the Intel TBB backend exposes a
# thread-count cap (tbb::global_control):
#   * libstdc++ (GCC, and Clang on most Linux distributions) always routes par through TBB.
#   * libc++ selects its PSTL backend at build time (TBB / std-thread / serial);
#     only the TBB build can be capped.
# CMake therefore locates or fetches TBB and exposes VASH_HAVE_TBB to the sources.
# The code itself distinguishes libstdc++ from libc++ via __GLIBCXX__ / _LIBCPP_VERSION.
#
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
#include <ciso646>
#if !defined(_LIBCPP_VERSION)
#error not libc++
#endif
int main() { return 0; }"
	VASH_STDLIB_IS_LIBCXX
)
if(VASH_STDLIB_IS_LIBCXX)
	message(STATUS "vash: C++ standard library is libc++")
else()
	message(STATUS "vash: C++ standard library is libstdc++ (par requires TBB)")
endif()

option(VASH_USE_TBB
	"Back parallel algorithms with Intel TBB and enable thread-count capping"
	ON
)
set(VASH_TBB_AVAILABLE FALSE)
if(VASH_USE_TBB)
	find_package(TBB CONFIG QUIET)
	if(TBB_FOUND)
		message(STATUS "vash: using system TBB ${TBB_VERSION}")
		set(VASH_TBB_AVAILABLE TRUE)
	else()
		message(STATUS "vash: system TBB not found, fetching oneTBB via FetchContent")
		# Build only the TBB runtime, not its own tests/examples.
		set(TBB_TEST OFF CACHE BOOL "" FORCE)
		set(TBB_STRICT OFF CACHE BOOL "" FORCE)
		set(TBBMALLOC_BUILD OFF CACHE BOOL "" FORCE)
		FetchContent_Declare(
			TBB
			GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB
			GIT_TAG        v2022.0.0   # adjust to a tag available upstream
			GIT_PROGRESS   TRUE
			SYSTEM
		)
		FetchContent_MakeAvailable(TBB)
		set(VASH_TBB_AVAILABLE TRUE)
	endif()
endif()
if(NOT VASH_TBB_AVAILABLE)
	message(WARNING
		"vash: building without TBB - parallel algorithms fall back to serial "
		"execution (std::execution::seq) and the thread-count ceiling is ignored."
	)
elseif(VASH_STDLIB_IS_LIBCXX)
	message(STATUS
		"vash: under libc++ the thread cap works only if its PSTL backend is TBB"
	)
endif()

include(GNUInstallDirs)
# executable
add_executable(ldblocks
	${CMAKE_CURRENT_LIST_DIR}/apps/ldblocks.cpp
)
target_include_directories(ldblocks
	PRIVATE ${brutilities_SOURCE_DIR}/include
	PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include
)
target_link_libraries(ldblocks
	PRIVATE vash
)
target_compile_options(ldblocks
	PRIVATE ${PROJECT_WARNINGS_CXX}
)
if(BUILD_TESTS AND SANITIZER_LIST)
	target_compile_options(ldblocks
		PRIVATE ${SANITIZER_COMPILE_FLAGS}
	)
	# Executables must also link the sanitizer runtime, not just compile with instrumentation.
	target_link_options(ldblocks
		PRIVATE -fsanitize=${SANITIZER_LIST}
	)
endif()
set_target_properties(ldblocks PROPERTIES
	CXX_STANDARD_REQUIRED ON
)
install(TARGETS ldblocks
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# library
add_library(vash
	${CMAKE_CURRENT_LIST_DIR}/src/gvarHash.cpp
	${CMAKE_CURRENT_LIST_DIR}/src/vashFunctions.cpp
	${CMAKE_CURRENT_LIST_DIR}/src/similarityMatrix.cpp
	${CMAKE_CURRENT_LIST_DIR}/src/vashLogging.cpp
	${brutilities_SOURCE_DIR}/src/random.cpp
)
target_include_directories(vash
	PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include
	PRIVATE ${brutilities_SOURCE_DIR}/include
)
set_target_properties(vash PROPERTIES
	POSITION_INDEPENDENT_CODE ON
	CXX_STANDARD_REQUIRED ON
)
target_compile_options(vash
	PRIVATE ${PROJECT_WARNINGS_CXX} -mbmi -mbmi2 -mpopcnt -lpthread
)
# Propagate the TBB backend (PUBLIC: consumers including the parallel shim need it too)
if(VASH_TBB_AVAILABLE)
	target_link_libraries(vash
		PUBLIC TBB::tbb
	)
	target_compile_definitions(vash
		PUBLIC VASH_HAVE_TBB
	)
endif()
if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTS AND SANITIZER_LIST)
	# vash is a static library (no link step), so it only needs the instrumentation flags;
	# the runtime is linked by the executables that consume it (ldblocks, tests).
	target_compile_options(vash
		PRIVATE ${SANITIZER_COMPILE_FLAGS}
	)
endif()
# Activate the in-source [vash-bench] phase timers when building benchmarks.
if(BUILD_BENCHMARKS)
	target_compile_definitions(vash
		PRIVATE VASH_BENCHMARK
	)
endif()
# no install for now

# tests
if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTS)
	add_executable(tests
		tests/tests.cpp
	)
	target_link_libraries(tests
		PRIVATE vash
		PRIVATE Catch2::Catch2WithMain
	)
	target_include_directories(tests
		PRIVATE include
		PRIVATE ${brutilities_SOURCE_DIR}/include
	)
	target_compile_options(tests
		PRIVATE -mpopcnt
	)
	if(SANITIZER_LIST)
		# Instrument the test translation unit and link the sanitizer runtime into the executable,
		# so the suite actually runs under the sanitizers (matching the instrumented vash library).
		target_compile_options(tests
			PRIVATE ${SANITIZER_COMPILE_FLAGS}
		)
		target_link_options(tests
			PRIVATE -fsanitize=${SANITIZER_LIST}
		)
	endif()
	include(CTest)
	include(Catch)
	catch_discover_tests(tests)
endif()

# benchmark harnesses (opt-in; activates VASH_BENCHMARK on the library above)
if(PROJECT_IS_TOP_LEVEL AND BUILD_BENCHMARKS)
	add_executable(ldGroupBench
		${CMAKE_CURRENT_LIST_DIR}/benchmarks/ldGroupBench.cpp
	)
	target_include_directories(ldGroupBench
		PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include
		PRIVATE ${brutilities_SOURCE_DIR}/include
	)
	target_link_libraries(ldGroupBench
		PRIVATE vash
	)
	target_compile_options(ldGroupBench
		PRIVATE ${PROJECT_WARNINGS_CXX}
	)
	set_target_properties(ldGroupBench PROPERTIES
		CXX_STANDARD_REQUIRED ON
	)
endif()
