cmake_minimum_required(VERSION 3.17)
project(SeniorDesign)

set(CMAKE_CXX_STANDARD 14)

add_executable(SeniorDesign CPU/main.cpp CPU/algorithm.cpp CPU/algorithm.h CPU/decoder.cpp
        CPU/decoder.h CPU/dependency.cpp CPU/dependency.h CPU/tree.cpp CPU/tree.h
        CPU/inputVerifier.cpp CPU/inputVerifier.h)