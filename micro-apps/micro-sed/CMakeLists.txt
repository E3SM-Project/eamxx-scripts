set(MICROSED_SRCS
  p3_ref.f90
  cpp_bridge.f90
  array_io.cpp
  initial_conditions.cpp
)

add_library(microsed ${MICROSED_SRCS})
target_include_directories (microsed PUBLIC ${SCREAM_MICRO_APPS_INCLUDES} ${SCREAM_MICRO_APPS_SOURCE_DIR}/share ${SCREAM_MICRO_APPS_SOURCE_DIR}/p3)

add_executable(unit_test_f90 unit_test.f90)
set_target_properties(unit_test_f90 PROPERTIES LINKER_LANGUAGE Fortran)
add_executable(run_and_cmp run_and_cmp.cpp)
add_executable(perf_cmp perf_cmp.cpp)

set(exes unit_test_f90 run_and_cmp perf_cmp)
set(p3_impls ref vanilla final)

foreach (p3_impl ${p3_impls})
  set(exe p3_${p3_impl})
  add_executable(${exe} ${exe}_driver.cpp)
  LIST(APPEND exes ${exe})
endforeach ()

foreach (exe ${exes})
  target_include_directories (${exe} PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
  target_link_libraries (${exe} share ${SCREAM_MICRO_APPS_LIBRARIES})
  target_link_libraries (${exe} microsed ${SCREAM_MICRO_APPS_LIBRARIES})
  if (NOT CUDA_BUILD)
    target_link_libraries (${exe} p3 ${SCREAM_MICRO_APPS_LIBRARIES})
  endif()
  set_target_properties (${exe} PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})
endforeach ()

add_custom_target(baseline
  COMMAND $<TARGET_FILE:run_and_cmp> -g run_and_cmp.baseline)

add_test(run_unit_test_f90 unit_test_f90)
add_test(p3_func_regression run_and_cmp run_and_cmp.baseline)
set(PERF_ARGS 10 128 300 30 1 3) # ni nk dt ts kdir repeat
foreach (p3_impl ${p3_impls})
  add_test(run_p3_${p3_impl} p3_${p3_impl} ${PERF_ARGS})
  if (NOT ${p3_impl} STREQUAL "ref")
    add_test(p3_perf_regr_${p3_impl} perf_cmp fortran_perf_run.dat ${p3_impl}_perf_run.dat)
  endif()
endforeach()
