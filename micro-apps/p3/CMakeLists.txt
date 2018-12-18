
# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  set(P3_SRCS
    p3_final.cpp
    p3_functions_find.cpp
    p3_functions_table3.cpp
    p3_functions_upwind.cpp
  )
  add_library(p3 ${P3_SRCS})
  target_include_directories (p3 PUBLIC ${SCREAM_MICRO_APPS_INCLUDES} ${SCREAM_MICRO_APPS_SOURCE_DIR}/share)
endif()

add_executable(p3_unit_test p3_unit_test.cpp)

set(exes p3_unit_test)

foreach (exe ${exes})
  target_include_directories (${exe} PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
  target_link_libraries (${exe} share ${SCREAM_MICRO_APPS_LIBRARIES})
  if (NOT CUDA_BUILD)
    target_link_libraries (${exe} p3 ${SCREAM_MICRO_APPS_LIBRARIES})
  endif()
  set_target_properties (${exe} PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})
endforeach ()

add_test(run_p3_unit_test p3_unit_test -u 16) # Cap at 16 threads
