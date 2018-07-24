set(P3_MICRO_SRCS
  micro_sed.f90
  array_io_mod.f90
  array_io.cpp
)

add_library(microsed ${P3_MICRO_SRCS})

add_executable(micro_sed micro_sed_driver.f90)
target_link_libraries(micro_sed microsed)

add_executable(unit_test_f90 unit_test.f90)
target_link_libraries(unit_test_f90 microsed)

add_executable(run_and_cmp run_and_cmp.cpp)
target_link_libraries(run_and_cmp microsed)

add_custom_target(baseline
  COMMAND $<TARGET_FILE:run_and_cmp> -g run_and_cmp.baseline)

add_executable(micro_sed_vanilla micro_sed_vanilla_driver.cpp)
target_link_libraries(micro_sed_vanilla microsed)

add_test(run_unit_test unit_test_f90)
add_test(run_micro_sed micro_sed 1 2 3 4 5 6 1)
add_test(micro_sed_func_regression run_and_cmp run_and_cmp.baseline)
add_test(run_micro_sed_vanilla micro_sed_vanilla 1 2 3 4 5 6 1)
