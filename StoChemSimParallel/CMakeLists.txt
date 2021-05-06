# Parallel CMakeLists

add_library(gpu_ssa STATIC
			cudaUtils.cpp
			cudaUtils.h 
			DirectMethod.cu 
			DirectMethod.cuh 
			inputParser.cpp 
			inputParser.h 
			kernels.cu 
			kernels.cuh 
			utils.cpp 
			utils.h 
			ArgParser.h
			ArgParser.cpp
			)

# Find boost and thrust (bundled with CUDA)
find_package(Boost 1.75.0 REQUIRED)
message("Found Boost include directories at ${Boost_INCLUDE_DIR}")
include_directories(${Boost_INCLUDE_DIR})

find_package(Thrust REQUIRED CONFIG)
thrust_create_target(Thrust)

# Set seperable compilation on the library and the executable
set_target_properties(gpu_ssa PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(gpu_ssa PROPERTIES CUDA_ARCHITECTURES 61)

# Tell CMake to compile Main with NVCC so that thrust can be compiled.
# set_source_files_properties(Main.cpp PROPERTIES LANGUAGE CUDA)
add_executable(stochemsim_gpu Main.cpp)

set_property(TARGET stochemsim_gpu PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET stochemsim_gpu PROPERTY CUDA_ARCHITECTURES 61)

# Link the boost, thrust (in CUDA), and in-house libraries to the executable

target_include_directories(stochemsim_gpu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

link_directories("${Boost_INCLUDE_DIR}/stage/lib/")
target_link_libraries(stochemsim_gpu gpu_ssa Thrust ${CUDA_curand_LIBRARY})
target_link_libraries(gpu_ssa Thrust)

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/banner_logo.txt"
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})