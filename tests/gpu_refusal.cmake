# A cook that refused for want of a GPU is a skip, not a failure, on a machine without one. ctest skips on the
# STARSCAPE-SKIP marker (SKIP_REGULAR_EXPRESSION), which only this file prints, so a refusal under
# STARSCAPE_REQUIRE_GPU (read when the test runs) fails instead: there, a device check that wrongly refuses every
# GPU must turn the gate red rather than skip it.
function(starscape_fail_run reason errors)
	string(FIND "${errors}" "no GPU to march on:" refused)
	if(NOT refused EQUAL -1 AND ("$ENV{STARSCAPE_REQUIRE_GPU}" STREQUAL "" OR "$ENV{STARSCAPE_REQUIRE_GPU}" STREQUAL "0"))
		message(FATAL_ERROR "STARSCAPE-SKIP (no GPU; set STARSCAPE_REQUIRE_GPU=1 to fail instead): ${errors}")
	endif()
	message(FATAL_ERROR "${reason}: ${errors}")
endfunction()
