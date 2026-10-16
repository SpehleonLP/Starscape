# Cooks SCENE with --threads 1 and --threads 8 and requires byte-identical sky and preview files.
# Pinned to the CPU: the GPU march has no threads, so under Auto this would test only the star draws.
foreach(threads 1 8)
	execute_process(COMMAND "${EXE}" cook "${SCENE}" -o "${OUT}/det${threads}.kre_starscape"
		--preview "${OUT}/det${threads}.hdr" --threads ${threads} --device cpu RESULT_VARIABLE code OUTPUT_QUIET)
	if(NOT code EQUAL 0)
		message(FATAL_ERROR "cook with ${threads} threads exited ${code}")
	endif()
endforeach()
foreach(extension kre_starscape hdr)
	execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${OUT}/det1.${extension}" "${OUT}/det8.${extension}" RESULT_VARIABLE differ)
	if(NOT differ EQUAL 0)
		message(FATAL_ERROR "${extension} differs between 1 and 8 threads")
	endif()
endforeach()
