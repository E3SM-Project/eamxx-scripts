add_executable (nanovec vec.cpp fvec.f90)
# Use a pack size that does not divide 128 to test correctness in that case.
target_compile_definitions (nanovec PRIVATE -DVEC_PACKN=13)
set (NANOVEC_DP 0)
if (${SCREAM_DOUBLE_PRECISION})
  set (NANOVEC_DP 1)
endif ()
target_compile_definitions (nanovec PRIVATE -DVEC_PACKN=13 -DVEC_DP=${NANOVEC_DP})
target_include_directories (nanovec PUBLIC ${SCREAM_MICRO_APPS_INCLUDES})
target_link_libraries (nanovec ${SCREAM_MICRO_APPS_LIBRARIES})
set_target_properties (nanovec PROPERTIES LINK_FLAGS ${SCREAM_MICRO_APPS_LINK_FLAGS})

add_test(nanovec nanovec -nc 4 -ns 100)
