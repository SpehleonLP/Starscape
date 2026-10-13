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
	src/population.cpp
	src/galaxy.cpp
	src/presets.cpp
	src/scene.cpp
	src/march.cpp
	src/cook.cpp
	src/writer.cpp
	src/preview.cpp
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
	tests/test_format.cpp
	tests/test_galaxy.cpp
	tests/test_population.cpp
	tests/test_scene.cpp
	tests/test_march.cpp
	tests/test_cook.cpp
	tests/test_writer.cpp
	tests/test_preview.cpp
)
add_executable(starscape-tests ${STARSCAPE_TESTS})
target_link_libraries(starscape-tests PRIVATE starscape-core)
add_test(NAME unit COMMAND starscape-tests)

add_executable(starscape src/cli.cpp)
target_link_libraries(starscape PRIVATE starscape-core)

set(CLI_OUT ${CMAKE_CURRENT_BINARY_DIR}/cli-test)
file(MAKE_DIRECTORY ${CLI_OUT})
file(WRITE ${CLI_OUT}/garbage.kre_starscape "this is not a sky")
# Extra checks go in the EXTRA list (-DSTDERR_CONTAINS=..., -DABSENT=...).
include(CMakeParseArguments)
function(starscape_cli_test name expect)
	cmake_parse_arguments(T "" "" "ARGV_LIST;EXTRA" ${ARGN})
	list(JOIN T_ARGV_LIST "|" joined)
	add_test(NAME ${name} COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:starscape> -DEXPECT=${expect} "-DARGS=${joined}" ${T_EXTRA} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
endfunction()
starscape_cli_test(cli-cook 0 ARGV_LIST cook ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/tiny.toml -o ${CLI_OUT}/tiny.kre_starscape --preview ${CLI_OUT}/tiny.hdr)
set_tests_properties(cli-cook PROPERTIES FIXTURES_SETUP tiny-sky)
starscape_cli_test(cli-inspect 0 ARGV_LIST inspect ${CLI_OUT}/tiny.kre_starscape)
starscape_cli_test(cli-preview 0 ARGV_LIST preview ${CLI_OUT}/tiny.kre_starscape -o ${CLI_OUT}/again.hdr --exposure 2)
set_tests_properties(cli-inspect cli-preview PROPERTIES FIXTURES_REQUIRED tiny-sky)
starscape_cli_test(cli-bad-scene 1 ARGV_LIST cook ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/bad-key.toml -o ${CLI_OUT}/bad.kre_starscape EXTRA -DSTDERR_CONTAINS=qualty -DABSENT=${CLI_OUT}/bad.kre_starscape)
starscape_cli_test(cli-missing-scene 2 ARGV_LIST cook ${CLI_OUT}/absent.toml -o ${CLI_OUT}/absent.kre_starscape)
starscape_cli_test(cli-garbage-file 1 ARGV_LIST inspect ${CLI_OUT}/garbage.kre_starscape EXTRA -DSTDERR_CONTAINS=garbage.kre_starscape)
starscape_cli_test(cli-missing-file 2 ARGV_LIST inspect ${CLI_OUT}/absent.kre_starscape)
starscape_cli_test(cli-usage 1 ARGV_LIST frobnicate)
starscape_cli_test(cli-draft 0 ARGV_LIST cook ${CMAKE_CURRENT_SOURCE_DIR}/scenes/milky-way-from-earth.toml -o ${CLI_OUT}/draft.kre_starscape --draft --preview ${CLI_OUT}/draft.hdr)
starscape_cli_test(cli-draft-stylized 0 ARGV_LIST cook ${CMAKE_CURRENT_SOURCE_DIR}/scenes/stylized-milky-way.toml -o ${CLI_OUT}/stylized.kre_starscape --draft --preview ${CLI_OUT}/stylized.hdr)

add_test(NAME cli-determinism COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:starscape> -DSCENE=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/tiny.toml
	-DOUT=${CLI_OUT} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_same.cmake)
