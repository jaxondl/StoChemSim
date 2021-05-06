#include "testfuncs.cuh"

__host__  void test_updateSimsAndPropensities() {
	const int s = 3;
	const int n = 3;
	const int m = 3;
	const int max_reactants = 2;
	const int max_affected = 3;
	double propensities[s * m] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double reaction_rates[n] = { 1, 2, 0.5 };
	// 2A -> B, 1
	// A 2B -> C, 2
	// C -> 2C, 0.5
	int reactants[n * max_reactants * 2] = { 0, 2, 0, 0, 0, 1, 1, 2, 2, 1, 0, 0 };
	int state_changes[m * n] = { -2, 1, 0, -1, -2, 1, 0, 0, 1 };
	int dep_matrix[n * max_affected] = { 0, 1, -1, 0, 1, 2, 2, -1, -1 };
	int sim_configs[m * s] = { 5, 5, 5, 5, 5, 5, 5, 5, 5 };
	int fired_reactions[s] = { 0, 1, 2 };
	bool stability_flags[s] = { false, false, false };

	double* dev_props = 0;
	double* dev_rates = 0;
	int* dev_reactants = 0;
	int* dev_state_changes = 0;
	int* dev_dep_matrix = 0;
	int* dev_fired_reactions = 0;
	int* dev_sim_configs = 0;
	bool* dev_stability_flags = 0;

	cudaMalloc(&dev_props, s * m * sizeof(double));
	cudaMalloc(&dev_rates, m * sizeof(double));
	cudaMalloc(&dev_reactants, n * max_reactants * 2 * sizeof(int));
	cudaMalloc(&dev_state_changes, m * n * sizeof(int));
	cudaMalloc(&dev_dep_matrix, n * max_affected * sizeof(int));
	cudaMalloc(&dev_sim_configs, m * s * sizeof(int));
	cudaMalloc(&dev_fired_reactions, s * sizeof(int));
	cudaMalloc(&dev_stability_flags, s * sizeof(bool));

	cudaMemcpy(dev_props, propensities, s * m * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_rates, reaction_rates, m * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_reactants, reactants, n * max_reactants * 2 * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_state_changes, state_changes, m * n * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_dep_matrix, dep_matrix, n * max_affected * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_sim_configs, sim_configs, m * s * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_fired_reactions, fired_reactions, s * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_stability_flags, stability_flags, s * sizeof(bool), cudaMemcpyHostToDevice);

	dim3 grid(1, 1);
	dim3 threads(16, 16);
	//updateSims(s, n, dev_fired_reactions, dev_sim_configs, dev_state_changes, dev_stability_flags);
	//updateProps(s, n, m, max_reactants, dev_sim_configs, dev_reactants, dev_rates, dev_props);

	cudaMemcpy(propensities, dev_props, s * m * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(reaction_rates, dev_rates, m * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(reactants, dev_reactants, n * max_reactants * 2 * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(state_changes, dev_state_changes, m * n * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(dep_matrix, dev_dep_matrix, n * max_affected * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(sim_configs, dev_sim_configs, m * s * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(fired_reactions, dev_fired_reactions, s * sizeof(int), cudaMemcpyDeviceToHost);

	std::cout << "Simulation Configs" << '\n';
	for (int i = 0; i < s * n; i++) {
		std::cout << sim_configs[i] << " ";
	}
	std::cout << "\n\n";

	std::cout << "Propensities" << '\n';
	for (int i = 0; i < s * m; i++) {
		std::cout << propensities[i] << " ";
	}
}

__host__ bool test_scaleRVs(int s, int m, bool verbose) {
	thrust::host_vector<double> propensities(s * m);
	for (int i = 0; i < s; i++) {
		for (int j = 0; j < m; j++) {
			propensities[i * m + j] = j + 1;
		}
	}

	// Construct keys for row-wise inclusive scan
	thrust::host_vector<int> keys(s * m);
	for (int i = 0; i < s; i++) {
		for (int j = 0; j < m; j++) {
			keys[i * m + j] = i; // key for each element in a row is just its row number
		}
	}

	thrust::host_vector<double> propscan(s * m);
	thrust::inclusive_scan_by_key(thrust::host, keys.begin(), keys.end(), propensities.begin(), propscan.begin());

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	thrust::host_vector<double> uniformRVs(s);
	for (int i = 0; i < s; i++) {
		uniformRVs[i] = distribution(generator);
	}

	thrust::host_vector<double> scaledRVs(s);
	for (int i = 0; i < s; i++) {
		scaledRVs[i] = uniformRVs[i] * propscan[m + i * m - 1];
	}

	thrust::device_vector<double> dev_propscan = propscan;
	thrust::device_vector<double> dev_uniformRVs = uniformRVs;
	//scaleRVs(thrust::raw_pointer_cast(dev_uniformRVs.data()), thrust::raw_pointer_cast(dev_propscan.data()), s, m);
	
	thrust::host_vector<double> par_scaledRVs = dev_uniformRVs;

	int fail_index = -1;
	for (int i = 0; i < s; i++) {
		if (par_scaledRVs[i] != scaledRVs[i]) {
			fail_index = i;
			break;
		}
	}

	if (fail_index != -1) {
		std::cout << "Failure for s = " << s << " and m = " << m << '\n';
		std::cout << "Fails first at index i = " << fail_index << "\n\n";
		if (verbose) {
			std::cout << "Random Variables" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << uniformRVs[i] << ", ";
			}
			std::cout << "}\n\n";

			std::cout << "Propensity Sums" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << propscan[m + i * m - 1] << ", ";
			}
			std::cout << "}\n\n";

			std::cout << "Host Result" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << scaledRVs[i] << ", ";
			}
			std::cout << "}\n\n";

			std::cout << "Device Result" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << par_scaledRVs[i] << ", ";
			}
			std::cout << "}\n\n";
		}
		return false;
	}
	else {
		std::cout << "Success for s = " << s << " and m = " << m << "\n\n";
		if (verbose) {
			std::cout << "Host Result" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << scaledRVs[i] << ", ";
			}
			std::cout << "}\n\n";

			std::cout << "Device Result" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << par_scaledRVs[i] << ", ";
			}
			std::cout << "}\n\n";
		}
		return true;
	}
}

__host__ bool test_checkBins(int s, int m, bool verbose) {
	thrust::host_vector<double> propensities(s * m);
	for (int i = 0; i < s; i++) {
		for (int j = 0; j < m; j++) {
			propensities[i * m + j] = j + 1;
		}
	}

	// Construct keys for row-wise inclusive scan
	thrust::host_vector<int> keys(s * m);
	for (int i = 0; i < s; i++) {
		for (int j = 0; j < m; j++) {
			keys[i * m + j] = i; // key for each element in a row is just its row number
		}
	}

	thrust::host_vector<double> propscan(s * m);
	thrust::inclusive_scan_by_key(thrust::host, keys.begin(), keys.end(), propensities.begin(), propscan.begin());

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	thrust::host_vector<double> uniformRVs(s);
	for (int i = 0; i < s; i++) {
		uniformRVs[i] = distribution(generator);
	}

	thrust::device_vector<double> dev_propscan = propscan;
	thrust::device_vector<double> dev_uniformRVs = uniformRVs;
	//scaleRVs(thrust::raw_pointer_cast(dev_uniformRVs.data()), thrust::raw_pointer_cast(dev_propscan.data()), s, m);

	// move scaled RVs back to host
	uniformRVs = dev_uniformRVs;

	thrust::host_vector<int> bins(s);
	thrust::device_vector<int> dev_bins(s);
	thrust::fill(dev_bins.begin(), dev_bins.end(), -1);

	// Compute bins on host
	for (int i = 0; i < s; i++) {
		double rv = uniformRVs[i];
		for (int j = 0; j < m; j++) {
			double left_edge, right_edge;
			if (j == 0) {
				left_edge = 0;
				right_edge = propscan[i * m + j];
			}
			else {
				left_edge = propscan[i * m + j - 1];
				right_edge = propscan[i * m + j];
			}

			// Push right edge over to include total for last element
			if (j == m - 1) {
				right_edge += 1;
			}

			if (left_edge <= rv && rv < right_edge) {
				bins[i] = j;
			};
		}
	}

	//checkBins(thrust::raw_pointer_cast(dev_propscan.data()), thrust::raw_pointer_cast(dev_uniformRVs.data()), thrust::raw_pointer_cast(dev_bins.data()), s, m);

	thrust::host_vector<int> par_bins = dev_bins;

	int fail_index = -1;
	for (int i = 0; i < s; i++) {
		if (par_bins[i] != bins[i]) {
			fail_index = i;
			break;
		}
	}

	if (fail_index != -1) {
		std::cout << "Failure for s = " << s << " and m = " << m << '\n';
		std::cout << "Fails first at index i = " << fail_index << "\n\n";
		if (verbose) {
			std::cout << "Random Variables" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << uniformRVs[i] << ", ";
			}
			std::cout << "}\n\n\n";

			std::cout << "Propensity Sums" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << propscan[m + i * m - 1] << ", ";
			}
			std::cout << "}\n\n";

			std::cout << "Host Result" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << bins[i] << ", ";
			}
			std::cout << "}\n\n";

			std::cout << "Device Result" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << par_bins[i] << ", ";
			}
			std::cout << "}\n\n";
		}
		return false;
	}
	else {
		std::cout << "Success for s = " << s << " and m = " << m << "\n\n";
		if (verbose) {
			std::cout << "Host Result" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << bins[i] << ", ";
			}
			std::cout << "}\n\n";

			std::cout << "Device Result" << '\n';
			std::cout << "{";
			for (int i = 0; i < s; i++) {
				std::cout << par_bins[i] << ", ";
			}
			std::cout << "}\n\n";
		}
		return true;
	}
}