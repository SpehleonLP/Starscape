cmake_minimum_required(VERSION 3.25)
project(Starscape VERSION 0.1 LANGUAGES C CXX)

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

# The GPU march: Vulkan loaded at run time through volk, so a binary without a Vulkan driver still runs
# (on the CPU); shaders compiled at build time and embedded, so it needs no files beside it.
find_program(GLSLANG_VALIDATOR glslangValidator)
if(GLSLANG_VALIDATOR)
	set(STARSCAPE_GPU_DEFAULT ON)
else()
	set(STARSCAPE_GPU_DEFAULT OFF)
endif()
option(STARSCAPE_GPU "Build the Vulkan compute march" ${STARSCAPE_GPU_DEFAULT})
if(STARSCAPE_GPU AND NOT GLSLANG_VALIDATOR)
	message(WARNING "STARSCAPE_GPU needs glslangValidator; building without the GPU march")
	set(STARSCAPE_GPU OFF)
endif()
if(NOT STARSCAPE_GPU)
	message(STATUS "Starscape: GPU march off; cooks run on the CPU")
endif()

set(SHADER_OUT ${CMAKE_CURRENT_BINARY_DIR}/shaders)
file(MAKE_DIRECTORY ${SHADER_OUT})
# Everything a shader may include, listed rather than globbed: a glob is read only when CMake configures, so a
# header added later would leave shaders built against a stale copy.
set(STARSCAPE_SHADER_INCLUDES
	src/field/common.glsl.h
	src/field/noise.glsl.h
	src/field/galaxy.glsl.h
	src/field/step.glsl.h
	shaders/face.glsl
	shaders/philox.glsl
)
list(TRANSFORM STARSCAPE_SHADER_INCLUDES PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)
# The embedded SPIR-V lands in <name>.spv.h: a plain <name>.h would shadow or be shadowed by src/<name>.h.
function(starscape_shader name)
	set(source ${CMAKE_CURRENT_SOURCE_DIR}/shaders/${name}.comp)
	set(field_headers ${STARSCAPE_SHADER_INCLUDES})
	add_custom_command(OUTPUT ${SHADER_OUT}/${name}.spv.h
		COMMAND ${GLSLANG_VALIDATOR} -V --target-env vulkan1.2 -I${CMAKE_CURRENT_SOURCE_DIR}/src -I${CMAKE_CURRENT_SOURCE_DIR}/shaders
			-o ${SHADER_OUT}/${name}.spv ${source}
		COMMAND ${CMAKE_COMMAND} -DIN=${SHADER_OUT}/${name}.spv -DOUT=${SHADER_OUT}/${name}.spv.h -DNAME=k_${name}_spirv
			-P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmbedSpirv.cmake
		DEPENDS ${source} ${field_headers} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmbedSpirv.cmake
		VERBATIM)
	set_property(GLOBAL APPEND PROPERTY STARSCAPE_SHADER_HEADERS ${SHADER_OUT}/${name}.spv.h)
endfunction()

if(STARSCAPE_GPU)
	add_library(starscape-volk STATIC third_party/volk/volk.c)
	target_include_directories(starscape-volk SYSTEM PUBLIC third_party/volk third_party/vulkan)
	target_compile_definitions(starscape-volk PUBLIC VK_NO_PROTOTYPES)
	target_link_libraries(starscape-volk PUBLIC ${CMAKE_DL_LIBS})
	target_sources(starscape-core PRIVATE src/gpu/device.cpp src/gpu/pack.cpp src/gpu/probe.cpp src/gpu/march.cpp)
	# Vulkan create-info structs are written with designated initializers that name only the fields in use;
	# the omitted ones are meant to be zero, which is what this warning objects to.
	if(NOT MSVC)
		set_source_files_properties(src/gpu/device.cpp src/gpu/probe.cpp src/gpu/march.cpp PROPERTIES COMPILE_OPTIONS -Wno-missing-field-initializers)
	endif()
	target_compile_definitions(starscape-core PUBLIC STARSCAPE_GPU=1)
	target_link_libraries(starscape-core PUBLIC starscape-volk)
	target_include_directories(starscape-core PUBLIC ${SHADER_OUT})
	starscape_shader(probe)
	starscape_shader(march)
	get_property(shader_headers GLOBAL PROPERTY STARSCAPE_SHADER_HEADERS)
	add_custom_target(starscape-shaders DEPENDS ${shader_headers})
	add_dependencies(starscape-core starscape-shaders)
else()
	target_sources(starscape-core PRIVATE src/gpu/march_stub.cpp)
endif()

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

if(STARSCAPE_GPU)
	add_executable(starscape-gpu-tests tests/gpu/main.cpp tests/gpu/test_device.cpp tests/gpu/test_probe.cpp tests/gpu/test_march.cpp)
	target_link_libraries(starscape-gpu-tests PRIVATE starscape-core)
	# test_scenes.h, shared with the unit tests.
	target_include_directories(starscape-gpu-tests PRIVATE tests)
	target_compile_definitions(starscape-gpu-tests PRIVATE STARSCAPE_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
	add_test(NAME gpu COMMAND starscape-gpu-tests)
	set_tests_properties(gpu PROPERTIES SKIP_RETURN_CODE 77)
endif()

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
starscape_cli_test(cli-device-bad 1 ARGV_LIST cook ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/tiny.toml -o ${CLI_OUT}/dev.kre_starscape --device tpu EXTRA -DSTDERR_CONTAINS=tpu)
# Each build's own refusal: the real device check's with the GPU march, the stub's without it.
if(STARSCAPE_GPU)
	set(gpu_refusal "the GPU is disabled by STARSCAPE_GPU_DISABLE")
else()
	set(gpu_refusal "this build has no GPU march")
endif()
starscape_cli_test(cli-device-gpu-disabled 1 ARGV_LIST cook ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/tiny.toml -o ${CLI_OUT}/devg.kre_starscape --device gpu EXTRA "-DSTDERR_CONTAINS=${gpu_refusal}")
starscape_cli_test(cli-device-auto-disabled 0 ARGV_LIST cook ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/tiny.toml -o ${CLI_OUT}/deva.kre_starscape --device auto EXTRA "-DSTDOUT_CONTAINS=marching on the CPU: ${gpu_refusal}")
set_tests_properties(cli-device-gpu-disabled cli-device-auto-disabled PROPERTIES ENVIRONMENT STARSCAPE_GPU_DISABLE=1)

add_test(NAME cli-determinism COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:starscape> -DSCENE=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/tiny.toml
	-DOUT=${CLI_OUT} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_same.cmake)

# A scene that sets no nebula key and one that writes each at its default must cook to the same bytes.
add_test(NAME cli-nebula-defaults-inert COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:starscape>
	-DSCENE_A=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nebula-absent.toml -DSCENE_B=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nebula-defaults.toml
	-DOUT=${CLI_OUT} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_same_cooks.cmake)
