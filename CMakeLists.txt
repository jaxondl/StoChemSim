# Root CMakeLists
cmake_minimum_required(VERSION 3.9)
project(StoChemSim LANGUAGES CXX)

option(TARGET_CUDA "Compile CUDA code" ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(sequential)
link_directories(${Boost_INCLUDE_DIR}/stage/lib/)

if(TARGET_CUDA)
	enable_language(CUDA)
	include_directories("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.2/include")
	add_subdirectory(parallel)
endif()
