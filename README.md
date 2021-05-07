# StoChemSim

Wolfram Language package to stochastically simulate chemical reaction networks.

Includes Gillespie's Direct SSA and Dr. Soloveichik's Bounded Tau Leaping algorithms.

Dependent on Dr. Soloveichik's CRNSimulator package.

## Wolfram Language Package Setup Instructions

**Note:** requires Wolfram Mathematica version 12.

**Installation instructions:** Drop StoChemSim.wl, CRNSimulator.m, and the "StoChemSimSequential" directory in the Mathematica Applications directory.
This will allow the package to be loaded using "<<StoChemSim`" in Mathematica.

The Mathematica Applications directory can be opened with "SystemOpen@FileNameJoin[{$UserBaseDirectory, "Applications"}]" in Mathematica.

Demo notebooks can be found in the WolframNotebooks directory.

## Sequential Command Line Tool Setup Instructions

**Step 1:** Ensure that TARGET_CUDA is set to OFF and then build with CMake.

**Step 2:** Run the executable. The command line arguments should be entered in the following order: executable name, input file path, output file path, tEnd, flag(s).

Example: ./stochemsim_direct ./inputs/large_crn.txt ./output.txt 0 -fo

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
