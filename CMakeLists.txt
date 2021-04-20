cmake_minimum_required(VERSION 3.13)
project(SeniorDesign)

set(CMAKE_CXX_STANDARD 14)

add_executable(SeniorDesign sequential/driverBTL.cpp sequential/directMethodSSA/directMethodSSA.cpp sequential/directMethodSSA/directMethodSSA.h sequential/common/decoder.cpp
        sequential/common/decoder.h sequential/directMethodSSA/dependencyGraph.cpp sequential/directMethodSSA/dependencyGraph.h sequential/directMethodSSA/reactionTree.cpp sequential/directMethodSSA/reactionTree.h
        sequential/common/inputVerifier.cpp sequential/common/inputVerifier.h parallel/common/gpuDecoder.cpp parallel/common/gpuDecoder.h sequential/boundedTauLeaping/boundedTauLeaping.h sequential/boundedTauLeaping/boundedTauLeaping.cpp)