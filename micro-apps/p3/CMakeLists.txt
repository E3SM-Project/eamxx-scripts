set(P3_MICRO_SRCS
  p3_ref.f90
  cpp_bridge.f90
  array_io.cpp
  initial_conditions.cpp
  util.cpp
  p3_final.cpp
  wsm.cpp
)

add_library(p3 ${P3_MICRO_SRCS})
target_include_directories (p3 PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})

add_executable(unit_test_f90 unit_test.f90)
set_target_properties(unit_test_f90 PROPERTIES LINKER_LANGUAGE Fortran)
add_executable(unit_test unit_test.cpp)
add_executable(run_and_cmp run_and_cmp.cpp)
add_executable(perf_cmp perf_cmp.cpp)

set(exes unit_test_f90 unit_test run_and_cmp perf_cmp)
set(p3_impls ref vanilla final)

foreach (p3_impl ${p3_impls})
  set(exe p3_${p3_impl})
  add_executable(${exe} ${exe}_driver.cpp)
  LIST(APPEND exes ${exe})
endforeach ()

foreach (exe ${exes})
  target_include_directories (${exe} PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
  target_link_libraries (${exe} p3 ${SCREAM_MICRO_APPS_LIBRARIES})
  set_target_properties (${exe} PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})
endforeach ()

add_executable (nanoapp vec.cpp fvec.f90)
# Use a pack size that does not divide 128 to test correctness in that case.
target_compile_definitions (nanoapp PRIVATE -DVEC_PACKN=13)
set (NANOAPP_DP 0)
if (${DOUBLE_PRECISION})
  set (NANOAPP_DP 1)
endif ()
target_compile_definitions (nanoapp PRIVATE -DVEC_PACKN=13 -DVEC_DP=${NANOAPP_DP})
target_include_directories (nanoapp PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
target_link_libraries (nanoapp ${SCREAM_MICRO_APPS_LIBRARIES})
set_target_properties (nanoapp PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})

add_custom_target(baseline
  COMMAND $<TARGET_FILE:run_and_cmp> -g run_and_cmp.baseline)

add_test(nanoapp nanoapp -nc 4 -ns 100)
add_test(run_unit_test_f90 unit_test_f90)
add_test(run_unit_test unit_test -u 16) # Cap at 16 threads
add_test(p3_func_regression run_and_cmp run_and_cmp.baseline)
set(PERF_ARGS 10 128 300 30 1 3) # ni nk dt ts kdir repeat
foreach (p3_impl ${p3_impls})
  add_test(run_p3_${p3_impl} p3_${p3_impl} ${PERF_ARGS})
  if (NOT ${p3_impl} STREQUAL "ref")
    add_test(p3_perf_regr_${p3_impl} perf_cmp fortran_perf_run.dat ${p3_impl}_perf_run.dat)
  endif()
endforeach()
