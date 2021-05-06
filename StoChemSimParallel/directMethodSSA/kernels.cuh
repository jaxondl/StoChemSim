#ifndef KERNELS_H
#define KERNELS_H

#include <cuda.h>
#include <curand.h>
#include <device_launch_parameters.h>

__global__ void updateTimesKernel(double* randomVariables, double* times, size_t s);

__global__ void offsetMultiplicationKernel(double* x1, double* x2, size_t n, size_t size);

__global__ void updateSimsKernel(int s, int n, int* fired_reactions, int* sim_configs, int* state_changes, bool* stability_flags);

__global__ void updatePropsKernel(int s, int n, int m, int max_reactants, int* sim_configs, int* reactants, double* reaction_rates, double* propensities);

__global__ void calculateExponentialRVsKernel(double* randomVariables, double* propscan, size_t s, size_t m);

__global__ void checkBinsKernel(double* propensity_scan, double* uniformRVs, int* bins, size_t s, size_t m);

__global__ void propCheckKernel(double* propscan, bool* stability_flags, int s, int m);

#endif