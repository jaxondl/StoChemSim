# StoChemSim

Wolfram Language package to stochastically simulate chemical reaction networks.

Includes Gillespie's Direct SSA and Dr. Soloveichik's Bounded Tau Leaping algorithms.

Dependent on Dr. Soloveichik's CRNSimulator package.

## Wolfram Language Package Setup Instructions

Installation instructions: Drop StoChemSim.wl, CRNSimulator.m, and the "StoChemSimSequential" directory in the Mathematica Applications directory.

This directory can be opened with "SystemOpen@FileNameJoin[{$UserBaseDirectory, "Applications"}]" in Mathematica.

Demo notebooks can be found in the WolframNotebooks directory.

## Sequential Command Line Tool Setup Instructions

**Step 1:** Using a command line front end, navigate to the StoChemSimSequential directory.

**Step 2:** Assuming C++ 11 is installed, compile the relevant .cpp files.

For direct SSA, enter the following command: g++ =std=c++11 driverSSA.cpp ./common/\*.cpp ./directMethodSSA/\*.cpp

For BTL, enter the following command: g++ =std=c++11 driverBTL.cpp ./common/\*.cpp ./boundedTauLeaping/\*.cpp

**Step 3:** Run the executable. The command line arguments should be entered in the following order: input file path, output file path, tEnd, flag(s).

Example: ./a.out ./inputs/large_crn.txt ./output.txt 0 -fo

## Parallel Command Line Tool Setup Instructions
**Step 1:** Download Boost version 1.75.0 and install it. Ensure that Boost's root is on your system's PATH variable.

**Step 2:** Download CUDA v11.2 and install it. Ensure that CUDA's root directory is given as CUDA_PATH in your system variables.

**Step 3:** Ensure that TARGET_CUDA is set to ON and then build with CMake.

**TROUBLESHOOTING** You may be missing CMake configuration files for Thrust and CUB. If CMake throws an error because of this, you can find the necessary files in StoChemSimParallel/external.

## Contributors

Ahad Ahmed

Tarek Allam

Isaac Lee

Jackson Lightfoot

Seth Sehon

Vidur Sinha

David Soloveichik

Zhecheng Wang
