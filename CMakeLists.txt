cmake_minimum_required(VERSION 3.25)
project(Starscape VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)

# Every source except the CLI entry point; the tests link the same library.
set(STARSCAPE_SOURCES
	src/noise.cpp
)

add_library(starscape-core STATIC ${STARSCAPE_SOURCES})
target_include_directories(starscape-core PUBLIC include src)
target_include_directories(starscape-core SYSTEM PUBLIC third_party)
target_compile_definitions(starscape-core PUBLIC TOML_EXCEPTIONS=0)
if(NOT MSVC)
	target_compile_options(starscape-core PUBLIC -Wall -Wextra -Wpedantic)
endif()
target_link_libraries(starscape-core PUBLIC Threads::Threads)

enable_testing()
set(STARSCAPE_TESTS
	tests/main.cpp
	tests/test_rng.cpp
	tests/test_noise.cpp
)
add_executable(starscape-tests ${STARSCAPE_TESTS})
target_link_libraries(starscape-tests PRIVATE starscape-core)
add_test(NAME unit COMMAND starscape-tests)
