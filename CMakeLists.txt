cmake_minimum_required(VERSION 3.23)
project(MatrixMultiplication)

set(CMAKE_CXX_STANDARD 14)

add_executable(MatrixMultiplication main.cpp ThreadedMatrixMult.h)
