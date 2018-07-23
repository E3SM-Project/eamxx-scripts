
set(P3_MICRO_SRCS
  micro_sed.f90
  array_io_mod.f90
  array_io.cpp
)

add_library(libmicrosed ${P3_MICRO_SRCS})

add_executable(micro_sed micro_sed_driver.f90)
target_link_libraries(micro_sed libmicrosed)

add_executable(unit_test_f90 unit_test.f90)
target_link_libraries(unit_test_f90 libmicrosed)

add_executable(run_and_cmp run_and_cmp.cpp)
target_link_libraries(run_and_cmp libmicrosed)

add_test(run_unit_test unit_test_f90)
add_test(run_micro_sed micro_sed 1 2 3 4 5 6 1)