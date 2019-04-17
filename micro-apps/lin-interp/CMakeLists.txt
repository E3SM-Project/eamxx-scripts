set(LI_SRCS
  li_ref.f90
  li_common.cpp
  li_cpp_bridge.f90
)

# LI CONFIG OPTIONS
set(LI_TIME_SETUP FALSE CACHE LOGICAL "Time the setup phase instead of LIs")

if (LI_TIME_SETUP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLI_TIME_SETUP")
endif()

add_library(li ${LI_SRCS})
target_include_directories (li PUBLIC ${SCREAM_MICRO_APPS_INCLUDES} ${SCREAM_MICRO_APPS_SOURCE_DIR}/share)

add_executable(li_unit_test li_unit_test.cpp)
add_executable(li_perf_cmp li_perf_cmp.cpp)

set(exes li_perf_cmp li_unit_test)
set(li_impls ref vanilla alg kokkos vect)

foreach (li_impl ${li_impls})
  set(exe li_${li_impl})
  add_executable(${exe} ${exe}_driver.cpp)
  LIST(APPEND exes ${exe})
endforeach ()

foreach (exe ${exes})
  target_include_directories (${exe} PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
  target_link_libraries (${exe} share ${SCREAM_MICRO_APPS_LIBRARIES})
  target_link_libraries (${exe} li ${SCREAM_MICRO_APPS_LIBRARIES})
#  if (NOT CUDA_BUILD)
#    target_link_libraries (${exe} li ${SCREAM_MICRO_APPS_LIBRARIES})
#  endif()
  set_target_properties (${exe} PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})
endforeach ()

# add_custom_target(baseline
#   COMMAND $<TARGET_FILE:run_and_cmp> -g run_and_cmp.baseline)

add_test(run_li_unit_test li_unit_test)
# add_test(p3_func_regression run_and_cmp run_and_cmp.baseline)
set(PERF_ARGS 10 128 128 0.01 10) # ncol km1 km2 minthresh repeat
foreach (li_impl ${li_impls})
   add_test(run_li_${li_impl} li_${li_impl} ${PERF_ARGS})
   if (NOT ${li_impl} STREQUAL "ref")
     add_test(li_perf_regr_${li_impl} li_perf_cmp fortran_perf_run.dat ${li_impl}_perf_run.dat)
   endif()
endforeach()
