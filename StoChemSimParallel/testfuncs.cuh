#ifndef TESTFUNCS_H
#define TESTFUNCS_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <cuda.h>
#include <device_launch_parameters.h>
#include <random>
#include <string>
#include <iostream>
#include "DirectMethod.cuh"

__host__ void test_updateSimsAndPropensities();

__host__ bool test_scaleRVs(int s, int m, bool verbose);

__host__ bool test_checkBins(int s, int m, bool verbose);

#endif