cmake_minimum_required(VERSION 3.15)
project(lab7)

set(CMAKE_CXX_STANDARD 98)

add_executable(lab7 main.cpp IMatOperations.h SimpleMatOperations.cpp SimpleMatOperations.h VectorizedMatOperations.cpp VectorizedMatOperations.h MatGenerator.cpp MatGenerator.h MatSize.h BLASMatOperations.cpp BLASMatOperations.h)