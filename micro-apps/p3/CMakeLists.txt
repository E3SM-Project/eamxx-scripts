set(P3_MICRO_SRCS
  micro_sed.f90
  cpp_bridge.f90
  array_io.cpp
  initial_conditions.cpp
  util.cpp
)

add_library(microsed ${P3_MICRO_SRCS})
target_include_directories (microsed PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})

add_executable(micro_sed micro_sed_driver.cpp)
target_link_libraries(micro_sed microsed)

add_executable(unit_test_f90 unit_test.f90)
target_link_libraries(unit_test_f90 microsed)
set_target_properties(unit_test_f90 PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(unit_test unit_test.cpp)
target_link_libraries(unit_test microsed)
target_include_directories (unit_test PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
target_link_libraries (unit_test ${SCREAM_MICRO_APPS_LIBRARIES})
set_target_properties (unit_test PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})

add_executable(run_and_cmp run_and_cmp.cpp)
target_link_libraries(run_and_cmp microsed)
target_include_directories (run_and_cmp PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
target_link_libraries (run_and_cmp ${SCREAM_MICRO_APPS_LIBRARIES})
set_target_properties (run_and_cmp PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})

add_custom_target(baseline
  COMMAND $<TARGET_FILE:run_and_cmp> -g run_and_cmp.baseline)

add_executable(micro_sed_vanilla micro_sed_vanilla_driver.cpp)
target_link_libraries(micro_sed_vanilla microsed)

add_executable(perf_cmp perf_cmp.cpp)
target_link_libraries(perf_cmp microsed)

foreach (lbl vanilla pack)
  add_executable(micro_sed_${lbl}_kokkos micro_sed_${lbl}_driver_kokkos.cpp)
  target_link_libraries(micro_sed_${lbl}_kokkos microsed)
  target_include_directories (micro_sed_${lbl}_kokkos PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
  target_link_libraries (micro_sed_${lbl}_kokkos ${SCREAM_MICRO_APPS_LIBRARIES})
  set_target_properties (micro_sed_${lbl}_kokkos PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})
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

add_executable(test_cuda test_cuda.cpp)
target_include_directories (test_cuda PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
target_link_libraries (test_cuda ${SCREAM_MICRO_APPS_LIBRARIES})
set_target_properties (test_cuda PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})

set(PERF_ARGS 10 128 300 30 1) # ni nk dt ts kdir
add_test(nanoapp nanoapp -nc 4 -ns 100)
add_test(run_unit_test_f90 unit_test_f90)
add_test(run_unit_test unit_test)
add_test(run_micro_sed micro_sed ${PERF_ARGS})
add_test(run_micro_sed_vanilla micro_sed_vanilla ${PERF_ARGS})
add_test(run_micro_sed_vanilla_kokkos micro_sed_vanilla_kokkos ${PERF_ARGS})
add_test(run_micro_sed_pack_kokkos micro_sed_pack_kokkos ${PERF_ARGS})
add_test(func_regression run_and_cmp run_and_cmp.baseline)
add_test(perf_regr_vanilla perf_cmp fortran_perf_run.dat vanilla_perf_run.dat)
add_test(perf_regr_kokkos perf_cmp fortran_perf_run.dat kokkos_perf_run.dat)
add_test(perf_regr_pack_kokkos perf_cmp fortran_perf_run.dat pack_kokkos_perf_run.dat)
