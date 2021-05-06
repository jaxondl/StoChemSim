# Root CMakeLists
cmake_minimum_required(VERSION 3.9)
project(StoChemSim LANGUAGES CXX)

option(TARGET_CUDA "Compile CUDA code" ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(StoChemSimSequential)
message("Found Boost .lib direcetory at ${Boost_INCLUDE_DIR}/stage/lib/")
link_directories("${Boost_INCLUDE_DIR}/stage/lib/")

if(TARGET_CUDA)
	enable_language(CUDA)
	message("Found CUDA includes at $ENV{CUDA_PATH}/include")
	include_directories("$ENV{CUDA_PATH}/include")
	add_subdirectory(StoChemSimParallel)
endif()
