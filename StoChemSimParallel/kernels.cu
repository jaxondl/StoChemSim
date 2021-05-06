#include "kernels.cuh"

// TESETER: Tarek
// Incremement times by drawn RVs
// double[] randomVariables: The random variables array on device
// double[] times: The array of times on device
// size_t s: the number of simulations
__global__ void updateTimesKernel(double* randomVariables, double* times, size_t s) {
	int gid = blockIdx.x * blockDim.x + threadIdx.x;
	if (gid < s) {
		times[gid] += randomVariables[gid];
	}
}

// TESTER: Tarek
// Multiply the elements of 1D array x1 by the last column of 2D array x2 in-place
// double[] x1: 1D Array which is multiplied and stores the output
// double[] x2: 2D Array which multiplies x1 in place
// int n: Offset for each row in order to grab the final entry
// int size: overall size of x1 and the number of rows in x2
__global__ void offsetMultiplicationKernel(double* x1, double* x2, size_t n, size_t size) {
	int gid = blockIdx.x * blockDim.x + threadIdx.x;

	// Ignore excess threads for a multi-block invocation
	if (gid < size) {
		x1[gid] *= x2[n + gid * n - 1];
	}
}

// TESTER: Vidur
// fired_reactions: list (length s) of the IDs of the reaction fired for each simulation
// sim_configs: matrix giving configuration for each simulation
// state_changes: State change vector for each reaction
__global__ void updateSimsKernel(int s, int n, int* fired_reactions, int* sim_configs, int* state_changes, bool* stability_flags) {
	// Add reactions vectors given by checkBins() to simulation configs
	int gid = blockIdx.y * n + blockIdx.x * blockDim.x + threadIdx.x;
	int rid = blockIdx.x * blockDim.x + threadIdx.x;

	// Only progress if within bounds and simulation has not reached stability
	if (rid < n && stability_flags[blockIdx.y] == false) {
		sim_configs[gid] = sim_configs[gid] + state_changes[(fired_reactions[blockIdx.y]) * n + rid];
	}
}

__global__ void updatePropsKernel(int s, int n, int m, int max_reactants, int* sim_configs, int* reactants, double* reaction_rates, double* propensities) {
	// Each block updates one propensity value
	int gid = blockIdx.y * m + blockIdx.x * blockDim.x + threadIdx.x;
	int rid = blockIdx.x * blockDim.x + threadIdx.x;

	if (rid < m) {
		int reaction_num = rid; // gid % m;
		int reactants_starting_idx = reaction_num * max_reactants * 2;
		double propensity = reaction_rates[reaction_num];
		for (int i = 0; i < max_reactants; i++) {
			int molecule_idx = i * 2 + reactants_starting_idx;
			int reactant_coef_idx = i * 2 + reactants_starting_idx + 1;
			int molecule_amt = sim_configs[blockIdx.y * n + reactants[molecule_idx]]; // [sim*n + reactants[molecule_idx]];
			for (int j = 0; j < reactants[reactant_coef_idx]; j++) {
				propensity *= (molecule_amt - j);
			}
		}
		propensities[gid] = propensity;
	}
}

// TESTER: Zhecheng
// helper function for exponential RVs
__global__ void calculateExponentialRVsKernel(double* randomVariables, double* propscan, size_t s, size_t m) {
	int gid = blockIdx.y * gridDim.x * blockDim.x + blockIdx.x * blockDim.x + threadIdx.x;
	if (gid < s) {
		randomVariables[gid] = -1 / propscan[m + gid * m - 1] * log(1 - randomVariables[gid]); // replace each uniform variable in first half of RV array with exponential counterpart
	}
}

// TESTER: Tarek
// Each row of blocks corresponds to a single simulation (gridDim.y = s)
// Each row has enough blocks to create enough threads to check every bin location at once (gridDim.x * blockDim.x >= m)
__global__ void checkBinsKernel(double* propensity_scan, double* uniformRVs, int* bins, size_t s, size_t m) {
	int tid = threadIdx.x;
	int rid = blockIdx.x * blockDim.x + tid;
	int gid = blockIdx.y * m + rid;

	// only threads with row index less than the number of bins will check so that excess threads in last block are excluded
	if (rid < m) {
		// Move uniform RV for simulation and the bins being checked to shared memory
		double urv = uniformRVs[blockIdx.y]; // make sure uniformRVs is pointer offset

		// TODO: Would be made more efficient by padding with 0 at the beginning of the scan, as this would eliminate the if statement
		// TODO: Each value is loaded twice, should be changed to avoid this
		double left_edge, right_edge;
		if (rid == 0) {
			left_edge = 0;
			right_edge = propensity_scan[gid];
		}
		else {
			left_edge = propensity_scan[gid - 1];
			right_edge = propensity_scan[gid];
		}

		// Last thread in row increases its right edge to include 1 in the boundary
		if (rid == m - 1) {
			right_edge += 1;
		}

		// Only one warp per row will diverge on this instruction
		if (left_edge <= urv && urv < right_edge) {
			bins[blockIdx.y] = rid;
		}
	}
}

__global__ void propCheckKernel(double* propscan, bool* stability_flags, int s, int m) {
	int tid = threadIdx.x;
	int gid = blockIdx.x * blockDim.x + tid;
	double prop_row_sum = propscan[m + tid * m - 1];

	// only set stability to true if the thread is within bounds and the simulation propensity sum is nonzero
	if (gid < s && prop_row_sum == 0) {
		stability_flags[gid] = true;
	}
}