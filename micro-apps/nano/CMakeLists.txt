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

add_test(nanoapp nanoapp -nc 4 -ns 100)
