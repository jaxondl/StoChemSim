#include "DirectMethod.cuh"

__host__ void calculateExponentialRVs(double* randomVariables, double* propscan, int s, int m, cudaStream_t stream) {
	int num_blocks = s / MAX_THREADS_PER_BLOCK;
	int remainder = s % MAX_THREADS_PER_BLOCK;

	if (num_blocks >= 1 && remainder == 0) {
		// if we have more elements than fit in a single block and the number of elements is evenly divisible by the number of threads per block, then we need only need one block per max_threads elements
		calculateExponentialRVsKernel <<<num_blocks, MAX_THREADS_PER_BLOCK, 0, stream >>> (randomVariables, propscan, s, m);
	}
	else if (num_blocks >= 1 && remainder != 0) {
		// if we have more elements than fit in a single block and the number of elements is not evenly divisible by the number of threads per block, we need an extra block to process the left overs
		calculateExponentialRVsKernel <<<num_blocks + 1, MAX_THREADS_PER_BLOCK, 0, stream >>> (randomVariables, propscan, s, m);
	}
	else {
		// Otherwise we only need one block because s < max threads so we initialize with s threads
		calculateExponentialRVsKernel <<<1, s, 0, stream >>> (randomVariables, propscan, s, m);
	}
}

// TESTER: Zhecheng
// the first half of the RV will be exponential
// the second half of the RVs will be uniform
__host__ void drawRVs(double* propscan, double* randomVariables, int s, int m, curandGenerator_t gen, cudaStream_t stream, bool states_only) {
	if (states_only) {
		curandGenerateUniformDouble(gen, randomVariables + s, s);
	}
	else {
		curandGenerateUniformDouble(gen, randomVariables, 2 * s); // params: curand RNG, ptr to output, size; range: [0, 1]
		calculateExponentialRVs(randomVariables, propscan, s, m, stream);
	}
}

// TESTER: Tarek
// double[] uniformRVs: array of uniform RVs on device
// double[] propensity_scan: scanned array of propensities on device
// int s: the number of simulations AKA the number of RVs drawn
// int m: the number of reactions AKA the offset used to multiply by the final column
__host__ void scaleRVs(double* uniformRVs, double* propensity_scan, int s, int m, cudaStream_t stream) {
	int num_blocks = s / MAX_THREADS_PER_BLOCK;
	int remainder = s % MAX_THREADS_PER_BLOCK;

	if (num_blocks >= 1 && remainder == 0) {
		// if we have more elements than fit in a single block and the number of elements is evenly divisible by the number of threads per block, then we need only need one block per max_threads elements
		offsetMultiplicationKernel<<<num_blocks, MAX_THREADS_PER_BLOCK, 0, stream >>>(uniformRVs, propensity_scan, m, s);
	}
	else if (num_blocks >= 1 && remainder != 0) {
		// if we have more elements than fit in a single block and the number of elements is not evenly divisible by the number of threads per block, we need an extra block to process the left overs
		offsetMultiplicationKernel<<<num_blocks + 1, MAX_THREADS_PER_BLOCK, 0, stream >>>(uniformRVs, propensity_scan, m, s);
	}
	else {
		// Otherwise we only need one block because s < max threads so we initialize with s threads
		offsetMultiplicationKernel<<<1, s, 0, stream >>>(uniformRVs, propensity_scan, m, s);
	}
}

// double[] propensity_scan: scanned propensities
// double[] uniformRVs: Uniform random variables scaled to match propensities
// int s: number of simulations
// int m: number of reactions
__host__ void checkBins(double* propensity_scan, double* uniformRVs, int* bins, int s, int m, cudaStream_t stream) {
	int num_blocks = m / MAX_THREADS_PER_BLOCK;
	int remainder = m % MAX_THREADS_PER_BLOCK;

	if (num_blocks >= 1 && remainder == 0) {
		// All threads in all blocks will be engaged, so no need for an extra block
		dim3 gsize(num_blocks, s);
		checkBinsKernel<<<gsize, MAX_THREADS_PER_BLOCK, 0, stream >>>(propensity_scan, uniformRVs, bins, s, m);
	}
	else if (num_blocks >= 1 && remainder != 0) {
		// Last block will not have all threads engaged, this is the extra block
		dim3 gsize(num_blocks + 1, s);
		checkBinsKernel<<<gsize, MAX_THREADS_PER_BLOCK, 0, stream >>>(propensity_scan, uniformRVs, bins, s, m);
	}
	else {
		// Entire calculation can be run in one block
		dim3 gsize(1, s);
		checkBinsKernel<<<gsize, m, 0, stream >>>(propensity_scan, uniformRVs, bins, s, m);
	}
}

__host__ void updateSims(int s, int n, int* fired_reactions, int* sim_configs, int* state_changes, bool* stability_flags, cudaStream_t stream) {
	int num_blocks = n / MAX_THREADS_PER_BLOCK;
	int remainder = n % MAX_THREADS_PER_BLOCK;

	if (num_blocks >= 1 && remainder == 0) {
		dim3 gsize(num_blocks, s);
		updateSimsKernel << <gsize, MAX_THREADS_PER_BLOCK, 0, stream >> > (s, n, fired_reactions, sim_configs, state_changes, stability_flags);
	}
	if (num_blocks >= 1 && remainder != 0) {
		dim3 gsize(num_blocks + 1, s);
		updateSimsKernel << <gsize, MAX_THREADS_PER_BLOCK, 0, stream >> > (s, n, fired_reactions, sim_configs, state_changes, stability_flags);
	}
	else {
		dim3 gsize(1, s);
		updateSimsKernel << <gsize, n, 0, stream >> > (s, n, fired_reactions, sim_configs, state_changes, stability_flags);
	}
}

__host__ void updateProps(int s, int n, int m, int max_reactants, int* sim_configs, int* reactants, double* reaction_rates, double* propensities, cudaStream_t stream) {
	int num_blocks = m / MAX_THREADS_PER_BLOCK;
	int remainder = m % MAX_THREADS_PER_BLOCK;

	if (num_blocks >= 1 && remainder == 0) {
		dim3 gsize(num_blocks, s);
		updatePropsKernel << <gsize, MAX_THREADS_PER_BLOCK, 0, stream >> > (s, n, m, max_reactants, sim_configs, reactants, reaction_rates, propensities);
	}
	if (num_blocks >= 1 && remainder != 0) {
		dim3 gsize(num_blocks + 1, s);
		updatePropsKernel << <gsize, MAX_THREADS_PER_BLOCK, 0, stream >> > (s, n, m, max_reactants, sim_configs, reactants, reaction_rates, propensities);
	}
	else {
		dim3 gsize(1, s);
		updatePropsKernel << <gsize, m, 0, stream >> > (s, n, m, max_reactants, sim_configs, reactants, reaction_rates, propensities);
	}
}

__host__ void updateTimes(double* randomVariables, double* times, int s, cudaStream_t stream) {
	int num_blocks = s / MAX_THREADS_PER_BLOCK;
	int remainder = s % MAX_THREADS_PER_BLOCK;

	if (num_blocks >= 1 && remainder == 0) {
		// if we have more elements than fit in a single block and the number of elements is evenly divisible by the number of threads per block, then we need only need one block per max_threads elements
		updateTimesKernel << <num_blocks, MAX_THREADS_PER_BLOCK, 0, stream >> > (randomVariables, times, s);
	}
	else if (num_blocks >= 1 && remainder != 0) {
		// if we have more elements than fit in a single block and the number of elements is not evenly divisible by the number of threads per block, we need an extra block to process the left overs
		updateTimesKernel << <num_blocks + 1, MAX_THREADS_PER_BLOCK, 0, stream >> > (randomVariables, times, s);
	}
	else {
		// Otherwise we only need one block because s < max threads so we initialize with s threads
		updateTimesKernel << <1, s, 0, stream >> > (randomVariables, times, s);
	}
}	

__host__ void propCheck(double* propscan, bool* stability_flags, int s, int m, cudaStream_t stream) {
	int num_blocks = s / MAX_THREADS_PER_BLOCK;
	int remainder = s % MAX_THREADS_PER_BLOCK;

	if (num_blocks >= 1 && remainder == 0) {
		propCheckKernel<<< num_blocks, MAX_THREADS_PER_BLOCK, 0, stream >> > (propscan, stability_flags, s, m);
	}
	if (num_blocks >= 1 && remainder != 0) {
		propCheckKernel<<< num_blocks + 1, MAX_THREADS_PER_BLOCK, 0, stream >> > (propscan, stability_flags, s, m);
	}
	else {
		propCheckKernel << <1, s, 0, stream>> > (propscan, stability_flags, s, m);
	}
}

__host__ thrust::host_vector<int> build_full_pci(std::vector<int> pci, int s) {
	thrust::host_vector<int> full_pci(s * pci.size());
	for (int idx = 0; idx < s; idx++) {
		thrust::copy(thrust::host, pci.begin(), pci.end(), &full_pci[idx * pci.size()]);
	}
	return full_pci;
}

// TODO: Multi-gpu support
// TODO: Make randomVariables resize when using --states-only

// Loads data into device memory, launches kernels, and runs until final iteration
// int[] state_change_matrix: flattened matrix whose rows are state change vectors corresponding to reactions by row-index and species by column-index
// double[] rrc_vector: an array of rate reaction constants corresponding to reactions by index
// int[] configuration_matrix: a flattened matrix of simulation configurations corresponding to simulations by row-index and species by column-index
// double[] propensity_matrix: a flattened matrix of propensities corresponding to simulations by row-index and reactions by column-index
// int[] reactants_table: a flattened reverse-lookup table; row-indices correspond to reactions, column-indices correspond to reactant id/count pairs; should be padded to the max number of reactants with (0, 0)
// int s: the number of simulations
// int n: the number of species
// int m: the number of reactions
// int max_reactants: the maximum number of reactants in any reaction
// int stop: the number of iterations to run
// bool verbose: Whether to print status messages to the console.
__host__ void directMethod(int* state_change_matrix, double* rrc_vector, int* configuration_matrix, double* propensity_matrix, 
	int* reactants_table, int s, int n, int m, int max_reactants, int stop, bool verbose, bool all_confs, bool stability_only, 
	bool early_stop, bool debug, bool states_only, std::string rng, unsigned long long seed, std::vector<int> pci) {

	bool partial = pci.size() > 0;
	std::string fext = ".bin";

	if (verbose) {
		std::cout << "Constructing variables and moving data to device..." << std::endl << std::flush;
	}

	// Construct keys for row-wise inclusive scan
	int* host_keys = new int[s * m];
	for (int i = 0; i < s; i++) {
		for (int j = 0; j < m; j++) {
			host_keys[i * m + j] = i; // key for each element in a row is just its row number
		}
	}

	// build full list of copy indices for --partial-copy
	thrust::host_vector<int> h_full_pci = build_full_pci(pci, s);

	// Move input data to device
	thrust::device_vector<int> scm(state_change_matrix, state_change_matrix + n * m);
	thrust::device_vector<double> rrc_vec(rrc_vector, rrc_vector + m);
	thrust::device_vector<int> confmat(configuration_matrix, configuration_matrix + s * n);
	thrust::device_vector<double> propmat(propensity_matrix, propensity_matrix + s * m);
	thrust::device_vector<int> reactants(reactants_table, reactants_table + m * max_reactants * 2);
	thrust::device_vector<int> keys(host_keys, host_keys + s * m);
	thrust::device_vector<double> times(s, 0); // all times are 0 to start
	
	// Allocate output locations on device
	thrust::device_vector<double> propscan(s * m);
	thrust::device_vector<double> randomVariables(2*s); // Contains both exponential and uniform RVs (in first and second halves respectively)
	thrust::device_vector<int> bins(s);
	thrust::device_vector<bool> stability_flags(s, false); // Each value indicates whether the simulation has reached stability (all propensities zero in that sim), assumed false at start
	thrust::device_vector<int> partial_confmat(s * pci.size());
	thrust::device_vector<int> full_pci = h_full_pci;

	// Allocate output locations on host with pinned memory
	pinnedBoolVector h_stability_flags(s, false);
	pinnedIntVector h_confmat(configuration_matrix, configuration_matrix + s * n);
	pinnedDoubleVector h_times(s, 0);
	
	// Fill partial confmat with initial partial config for record keeping
	pinnedIntVector h_partial_confmat(s * pci.size());
	if (partial) {
		thrust::gather(thrust::host, h_full_pci.begin(), h_full_pci.end(), configuration_matrix, h_partial_confmat.begin());
	}

	// Pinned memory pointers for saving binaries and early stopping
	bool* h_stability_flags_ptr = thrust::raw_pointer_cast(h_stability_flags.data());
	int* h_confmat_ptr = thrust::raw_pointer_cast(h_confmat.data());
	double* h_times_ptr = thrust::raw_pointer_cast(h_times.data());
	int* h_partial_confmat_ptr = thrust::raw_pointer_cast(h_partial_confmat.data());

	// Raw device pointers for kernel processing outside thrust
	int* scm_ptr = thrust::raw_pointer_cast(scm.data());
	int* confmat_ptr = thrust::raw_pointer_cast(confmat.data());
	int* reactants_ptr = thrust::raw_pointer_cast(reactants.data());
	int* bins_ptr = thrust::raw_pointer_cast(bins.data());
	int* partial_confmat_ptr = thrust::raw_pointer_cast(partial_confmat.data());

	double* rrc_ptr = thrust::raw_pointer_cast(rrc_vec.data());
	double* uniformRVs_ptr = thrust::raw_pointer_cast(randomVariables.data()) + s;
	double* exponentialRVs_ptr = thrust::raw_pointer_cast(randomVariables.data());
	double* propmat_ptr = thrust::raw_pointer_cast(propmat.data());
	double* propscan_ptr = thrust::raw_pointer_cast(propscan.data());
	double* times_ptr = thrust::raw_pointer_cast(times.data());

	bool* stability_flags_ptr = thrust::raw_pointer_cast(stability_flags.data());

	if (verbose) {
		std::cout << "Setting stream policy..." << std::endl << std::flush;
	}

	// Create streams for kernel launches and memcpyasync (default stream used if not set)
	cudaStream_t kernel_stream, memcpy_stream;
	if (all_confs || stability_only || early_stop) {
		cudaStreamCreate(&kernel_stream);
		cudaStreamCreate(&memcpy_stream);
	}
	else {
		kernel_stream = 0; // 0 is default stream
		memcpy_stream = 0;
	}

	if (verbose) {
		std::cout << "Building random number generator..." << std::endl << std::flush;
	}

	// TODO: Fix exponential random variables
	// Make curand generator
	curandGenerator_t gen;
	curandRngType_t gen_type;
	
	if (rng == "XORWOW") {
		gen_type = CURAND_RNG_PSEUDO_XORWOW;
	}
	else if (rng == "MRG32K3A") {
		gen_type = CURAND_RNG_PSEUDO_MRG32K3A;
	}
	else if (rng == "MTGP32") {
		gen_type = CURAND_RNG_PSEUDO_MTGP32;
	}
	else if (rng == "PHILOX-4X32-10") {
		gen_type = CURAND_RNG_PSEUDO_PHILOX4_32_10;
	}
	else if (rng == "SOBOL32") {
		gen_type = CURAND_RNG_QUASI_SOBOL32;
	}
	else if (rng == "SOBOL64") {
		gen_type = CURAND_RNG_QUASI_SOBOL64;
	}
	else if (rng == "scrambledSOBOL32") {
		gen_type = CURAND_RNG_QUASI_SCRAMBLED_SOBOL32;
	}
	else if (rng == "scrambledSOBOL64") {
		gen_type = CURAND_RNG_QUASI_SCRAMBLED_SOBOL64;
	}
	else {
		gen_type = CURAND_RNG_PSEUDO_MT19937;
	}
	
	curandCreateGenerator(&gen, gen_type); // create MT19937 cuRAND RNG
	curandSetPseudoRandomGeneratorSeed(gen, seed);	// gen, seed

	// Make curand execute on kernel stream
	curandSetStream(gen, kernel_stream);

	// Make timestamp to ensure unique output binary filenames
	time_t rawtime;
	struct tm* timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, sizeof(buffer), "%d-%m-%Y_%H-%M-%S", timeinfo);
	std::string timestamp(buffer);

	// Save initial config for consistency, prior to main loop since it is not stored yet
	if (all_confs) {
		if (partial) {
			save_config("out/config0_" + timestamp + fext, std::vector<int>(h_partial_confmat.begin(), h_partial_confmat.end()));
		}
		else {
			save_config("out/config0_" + timestamp + fext, std::vector<int>(state_change_matrix, state_change_matrix + s * n));
		}

		if (!states_only) {
			save_times("out/times0_" + timestamp + fext, std::vector<double>(s * n, 0));
		}
	}

	// Only used if debug is set
	thrust::host_vector<double> h_propscan(s * m);
	thrust::host_vector<double> h_props(s * m);
	thrust::host_vector<bool> h_stabflags(s);
	thrust::host_vector<int> h_bins(s);
	thrust::host_vector<double> h_rvs(2*s);
	thrust::host_vector<double> h_chrono(s);
	thrust::host_vector<int> h_gather_out(s * pci.size());
	thrust::host_vector<int> h_configmat(s * n);

	// if not debugging deallocate to save space.
	if (!debug) {
		h_propscan.clear();
		h_propscan.shrink_to_fit();

		h_props.clear();
		h_props.shrink_to_fit();

		h_stabflags.clear();
		h_stabflags.shrink_to_fit();

		h_bins.clear();
		h_bins.shrink_to_fit();

		h_rvs.clear();
		h_rvs.shrink_to_fit();

		h_chrono.clear();
		h_chrono.shrink_to_fit();
	}

	int idx = 0;
	while(idx < stop || stability_only) {
		if (verbose && !debug) {
			if (stability_only) {
				std::cout << "\rRunning iteration [" << idx + 1 << "]" << std::flush;
			}
			else {
				std::cout << "\rRunning iteration [" << idx + 1 << "/" << stop << "]" << std::flush;
			}
		}

		// Update propensities
		updateProps(s, n, m, max_reactants, confmat_ptr, reactants_ptr, rrc_ptr, propmat_ptr, kernel_stream);

		// Inclusive scan of each row
		thrust::inclusive_scan_by_key(thrust::cuda::par.on(kernel_stream), keys.begin(), keys.end(), propmat.begin(), propscan.begin());

		// Set stability flags
		propCheck(propscan_ptr, stability_flags_ptr, s, m, kernel_stream);

		// Draw RVs
		drawRVs(propscan_ptr, exponentialRVs_ptr, s, m, gen, kernel_stream, states_only);

		// Update times
		updateTimes(exponentialRVs_ptr, times_ptr, s, kernel_stream); //TODO: stop time update once stability has been reached

		// Scale uniform RVs with element-wise multiplication
		scaleRVs(uniformRVs_ptr, propscan_ptr, s, m, kernel_stream);

		// Checkbins kernel
		checkBins(propscan_ptr, uniformRVs_ptr, bins_ptr, s, m, kernel_stream);

		if (debug) {
			// test output
			h_propscan = propscan;
			h_props = propmat;
			h_stabflags = stability_flags;
			h_bins = bins;
			h_rvs = randomVariables;
			h_chrono = times;
			h_gather_out = h_partial_confmat;
			h_configmat = confmat;
		}

		// Update state
		updateSims (s, n, bins_ptr, confmat_ptr, scm_ptr, stability_flags_ptr, kernel_stream);

		// Gather desired counts into temporary buffer
		if (partial) {
			thrust::gather(thrust::cuda::par.on(kernel_stream), full_pci.begin(), full_pci.end(), confmat.begin(), partial_confmat.begin());
		}

		// Let kernels get loaded before running blocking host code

		// Synchronize with device so that memory transfers from last round can finish executing.
		cudaDeviceSynchronize();

		if (debug) {
			std::cout << std::endl << std::endl;
			std::cout << "Iteration " << idx << std::endl;

			std::cout << "Stability Flags = {";
			for (int i = 0; i < s; i++) {
				std::cout << h_stabflags[i] << ", ";
			}
			std::cout << "}" << std::endl << std::endl;

			std::cout << "Configuration Matrix = {";
			for (int i = 0; i < s; i++) {
				std::cout << "{";
				for (int j = 0; j < n; j++) {
					std::cout << h_configmat[i * n + j] << ", ";
				}
				std::cout << "}, " << std::endl;
			}
			std::cout << "}" << std::endl;

			std::cout << "Exponential RVs = {";
			for (int i = 0; i < s; i++) {
				std::cout << h_rvs[i] << ", ";
			}
			std::cout << "}" << std::endl;

			std::cout << "Time Variables = {";
			for (int i = 0; i < s; i++) {
				std::cout << h_chrono[i] << ", ";
			}
			std::cout << "}" << std::endl << std::endl;

			std::cout << "Propensity Sums = {";
			for (int i = 0; i < s; i++) {
				std::cout << h_propscan[m + i * m - 1] << ", ";
			}
			std::cout << "}" << std::endl;

			std::cout << "Uniform RVs = {";
			for (int i = 0; i < s; i++) {
				std::cout << h_rvs[s + i] << ", ";
			}
			std::cout << "}" << std::endl;

			std::cout << "Propensity bounds = {";
			for (int i = 0; i < s; i++) {
				std::cout << "(" << h_propscan[i*m + h_bins[i] - 1] << ", " << h_propscan[i*m + h_bins[i]] << "), ";
			}
			std::cout << "}" << std::endl << std::endl;

			if (partial) {
				std::cout << "Partial config = {";
				for (int i = 0; i < pci.size(); i++) {
					std::cout << h_gather_out[i] << ", ";
				}
				std::cout << "}" << std::endl;
			}
		}
		
		// Don't run saving and checks on the first iteration since no device data hasn't been copied back yet
		if (idx > 0) {
			if (all_confs) {
				std::string fname = "out/config" + std::to_string(idx) + "_" + timestamp + fext;
				if (partial) {
					h_partial_confmat = partial_confmat;
					save_config(fname, std::vector<int>(h_partial_confmat_ptr, h_partial_confmat_ptr + s * pci.size()));
				}
				else {
					h_confmat = confmat;
					save_config(fname, std::vector<int>(h_confmat_ptr, h_confmat_ptr + s * n));
				}
				fname = "out/times" + std::to_string(idx) + "_" + timestamp + fext;

				if (!states_only) {
					save_times(fname, std::vector<double>(h_times_ptr, h_times_ptr + s));
				}
			}
			bool stability_reached = is_stable(h_stability_flags_ptr, s);
			if ((stability_only || early_stop) && stability_reached) {
				break;
			}
		}

		// asynchronously move data from device to host in seperate stream from kernels.
		// Should only not occur if we are running in vanilla, iterations-only, no-saving mode.
		if (all_confs || early_stop || stability_only) {
			if (partial) {
				cudaMemcpyAsync(h_partial_confmat_ptr, partial_confmat_ptr, s * pci.size() * sizeof(int), cudaMemcpyDeviceToHost, memcpy_stream);
			}
			else {
				cudaMemcpyAsync(h_confmat_ptr, confmat_ptr, s * n * sizeof(int), cudaMemcpyDeviceToHost, memcpy_stream);
			}
			cudaMemcpyAsync(h_stability_flags_ptr, stability_flags_ptr, s * sizeof(bool), cudaMemcpyDeviceToHost, memcpy_stream);
		}

		++idx;
	}

	if (!debug) {
		std::cout << std::endl << std::endl;
	}

	// Allow final memmory transfer to occur before continuing, if it is still running.
	cudaDeviceSynchronize();

	// Transfer of final config and save it
	if (partial) {
		h_partial_confmat = partial_confmat;
		save_config("out/config" + std::to_string(idx) + "_" + timestamp + fext, std::vector<int>(h_partial_confmat_ptr, h_partial_confmat_ptr + s * pci.size()));
	}
	else {
		h_confmat = confmat;
		save_config("out/config" + std::to_string(idx) + "_" + timestamp + fext, std::vector<int>(h_confmat_ptr, h_confmat_ptr + s * n));
	}

	if (!states_only) {
		h_times = times;
		save_times("out/times" + std::to_string(idx) + "_" + timestamp + fext, std::vector<double>(h_times_ptr, h_times_ptr + s));
	}

	if (verbose) {
		std::cout << "Simulation complete." << std::endl;
	}
}