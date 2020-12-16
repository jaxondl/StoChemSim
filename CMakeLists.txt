cmake_minimum_required(VERSION 3.13)
project(SeniorDesign)

set(CMAKE_CXX_STANDARD 14)

add_executable(SeniorDesign CPU/main.cpp CPU/directMethodSSA.cpp CPU/directMethodSSA.h CPU/decoder.cpp
        CPU/decoder.h CPU/dependencyGraph.cpp CPU/dependencyGraph.h CPU/reactionTree.cpp CPU/reactionTree.h
        CPU/inputVerifier.cpp CPU/inputVerifier.h)