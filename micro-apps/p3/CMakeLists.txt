
set(P3_MICRO_SRCS
  array_io_mod.f90
  array_io.cpp
)

add_executable(micro_sed micro_sed.f90 ${P3_MICRO_SRCS})
add_executable(unit_test_f90 unit_test.f90 ${P3_MICRO_SRCS})

add_test(run_unit_test unit_test_f90)
add_test(run_micro_sed micro_sed)
