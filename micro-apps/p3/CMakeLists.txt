
set(P3_MICRO_SRCS
  micro_sed.f90
)

add_executable(micro_sed ${P3_MICRO_SRCS})

add_test(run_micro_sed micro_sed)
