# Cooks SCENE in two separate runs, with --threads 1 and --threads 8, and requires byte-identical sky and preview
# files: the GPU march must not vary between runs, and the threaded calibration and star draws must not vary with the
# thread count. Earlier outputs are removed first, so a failed cook can never pass on stale files.
foreach(threads 1 8)
	file(REMOVE "${OUT}/det${threads}.kre_starscape" "${OUT}/det${threads}.hdr")
	execute_process(COMMAND "${EXE}" cook "${SCENE}" -o "${OUT}/det${threads}.kre_starscape"
		--preview "${OUT}/det${threads}.hdr" --threads ${threads} RESULT_VARIABLE code OUTPUT_QUIET ERROR_VARIABLE errors)
	if(NOT code EQUAL 0)
		message(FATAL_ERROR "cook with ${threads} threads exited ${code}: ${errors}")
	endif()
endforeach()
foreach(extension kre_starscape hdr)
	execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${OUT}/det1.${extension}" "${OUT}/det8.${extension}" RESULT_VARIABLE differ)
	if(NOT differ EQUAL 0)
		message(FATAL_ERROR "${extension} differs between two runs (1 and 8 threads)")
	endif()
endforeach()
