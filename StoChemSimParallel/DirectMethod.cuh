#ifndef DIRECT_METHOD_H
#define DIRECT_METHOD_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/gather.h>
#include <cuda.h>
#include <device_launch_parameters.h>
#include <curand.h>
#include <ctime>
#include <iostream>

///test 
#include <numeric>

#include "cudaUtils.h"
#include "kernels.cuh"

#define MAX_THREADS_PER_BLOCK 512

// Pinned allocators for thrust
typedef double mydouble;
typedef int myint;
typedef bool mybool;

typedef thrust::host_vector<mydouble, thrust::cuda::experimental::pinned_allocator<mydouble>> pinnedDoubleVector;
typedef thrust::host_vector<myint, thrust::cuda::experimental::pinned_allocator<myint>> pinnedIntVector;
typedef thrust::host_vector<mybool, thrust::cuda::experimental::pinned_allocator<mybool>> pinnedBoolVector;


__host__ void calculateExponentialRVs(double* randomVariables, double* propscan, int s, int m, cudaStream_t stream);

__host__ void drawRVs(double* propscan, double* randomVariables, int s, int m, curandGenerator_t gen, cudaStream_t stream, bool states_only);

__host__ void scaleRVs(double* uniformRVs, double* propensity_scan, int s, int m, cudaStream_t stream);

__host__ void checkBins(double* propensity_scan, double* uniformRVs, int* bins, int s, int m, cudaStream_t stream);

__host__ void updateSims(int s, int n, int* fired_reactions, int* sim_configs, int* state_changes, bool* stability_flags, cudaStream_t stream);

__host__ void updateProps(int s, int n, int m, int max_reactants, int* sim_configs, int* reactants, double* reaction_rates, double* propensities, cudaStream_t stream);

__host__ void updateTimes(double* randomVariables, double* times, int s, cudaStream_t stream);

__host__ void propCheck(double* propscan, bool* stability_flags, int s, int m, cudaStream_t stream);

__host__ thrust::host_vector<int> build_full_pci(thrust::host_vector<int> pci, int s);

__host__ void directMethod(int* state_change_matrix, double* rrc_vector, int* configuration_matrix, double* propensity_matrix, int* reactants_table, 
	int s, int n, int m, int max_reactants, int stop, bool verbose, bool all_confs, bool stability_only, bool early_stop, bool debug, bool states_only,
	std::string rng, unsigned long long seed, std::vector<int> pci);

#endif