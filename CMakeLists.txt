cmake_minimum_required(VERSION 3.13)
project(SeniorDesign)

set(CMAKE_CXX_STANDARD 14)

add_executable(SeniorDesign StoChemSimSequential/driverBTL.cpp StoChemSimSequential/directMethodSSA/directMethodSSA.cpp StoChemSimSequential/directMethodSSA/directMethodSSA.h StoChemSimSequential/common/decoder.cpp
        StoChemSimSequential/common/decoder.h StoChemSimSequential/directMethodSSA/dependencyGraph.cpp StoChemSimSequential/directMethodSSA/dependencyGraph.h StoChemSimSequential/directMethodSSA/reactionTree.cpp StoChemSimSequential/directMethodSSA/reactionTree.h
        StoChemSimSequential/common/inputVerifier.cpp StoChemSimSequential/common/inputVerifier.h StoChemSimParallel/common/gpuDecoder.cpp StoChemSimParallel/common/gpuDecoder.h StoChemSimSequential/boundedTauLeaping/boundedTauLeaping.h StoChemSimSequential/boundedTauLeaping/boundedTauLeaping.cpp)