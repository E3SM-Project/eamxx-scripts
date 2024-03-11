#!/usr/bin/bash

SRC_DIR=./

rm -rf CMake*

cmake \
 -D CMAKE_BUILD_TYPE:STIRNG=DEBUG     \
 -D CMAKE_CXX_COMPILER:STRING=mpicxx  \
 -D CMAKE_CXX_COMPILER:STRING=mpicxx  \
 ${SRC_DIR}
