cmake_minimum_required(VERSION 3.11)
project(vash VERSION 0.1 LANGUAGES CXX)

# Adding a Profile build type
get_property(isMultiConfig GLOBAL
	PROPERTY GENERATOR_IS_MULTICONFIG
)
if(isMultiConfig)
	if(NOT "Profile" IN_LIST CMAKE_CONFIGURATION_TYPES)
		list(APPEND CMAKE_CONFIGURATION_TYPES Profile)
	endif()
else()
	set(allowedBuildTypes Debug Release Profile)
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
set(CMAKE_CXX_FLAGS_PROFILE
	"-g -p -O3"
)
set(CMAKE_CXX_SHARED_LINKER_FLAGS_PROFILE
	"-g -p -O3"
)
set(CMAKE_CXX_STATIC_LINKER_FLAGS_PROFILE
	"-g -p -O3"
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
	PRIVATE -Wall -Wextra -pedantic
)
set_target_properties(ldblocks PROPERTIES
	CXX_STANDARD         14
	CXX_STANDARD_REQURED ON
	CXX_EXTENSIONS       OFF
)
install(TARGETS ldblocks
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# library
add_library(vash
	src/gvarHash.cpp
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
	PRIVATE -Wall -Wextra -pedantic -pthread
)
set_target_properties(vash PROPERTIES
	CXX_STANDARD         14
	CXX_STANDARD_REQURED ON
	CXX_EXTENSIONS       OFF
)
# no install for now

