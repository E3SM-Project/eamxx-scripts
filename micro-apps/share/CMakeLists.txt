set(SHARE_SRCS
  util.cpp
  scream_arch.cpp
)

add_library(share ${SHARE_SRCS})
target_include_directories (share PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})

add_executable(share_unit_test share_unit_test.cpp)

set(exes share_unit_test)

foreach (exe ${exes})
  target_include_directories (${exe} PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
  target_link_libraries (${exe} share ${SCREAM_MICRO_APPS_LIBRARIES})
  set_target_properties (${exe} PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})
endforeach ()

add_test(run_share_unit_test share_unit_test -u 16) # Cap at 16 threads
