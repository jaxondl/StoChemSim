#pragma warning( disable : 4267)

#include "inputParser.h"
#include "DirectMethod.cuh"
#include "ArgParser.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

#define MAX_SIMS 65535 // This is the maximum dimension of the y/z dimension of a cuda block grid (2^31 - 1). TODO: Increase by using the z dimension for simulations as well.

// Program options
// -h --help = output help string
// -f --infile = address of the input file being processed
// -s --sims = integer number of simulations to run
// -i --iters = number of iterations to run at most
// -v --verbose = run with logging to stdout (default false)
// -t --time = whether or not to record execution time information with --verbose (default false)
// -a --all-confs = return all configurations (default false)
// -o --stability-only = stop only once stability is reached if true (default false)
// -e --early-stop = stop when stability is reached if true (default false)
// -r --random-seed = integer seed to use for RNG (default 1234)
// --states-only = whether or not to record just the states or times as well (default false)
// --rng = one of several strings (see help), indicates which rng to use (default 'MT19937') 
// --partial-copy = a list of integer indicating the species indices to copy, rather than the whole confmat
// -g --gpus = number of GPUs to use (default 1) [NI]
int main(int argc, char* argv[]) {
	ArgParser* parser = new ArgParser(argc, argv);

	arg_status_t status = parser->parse();

	if (status == ERROR || status == HELP) {
		std::cout << std::endl;
		return 0;
	}

	std::string infile = parser->getInfile();
	std::string rng = parser->getRNG();
	std::vector<int> partial_copy_indices = parser->getPartialCopyIndices();
	int sims = parser->getSims();
	int iters = parser->getIters();
	int max_reactants = parser->getMaxReactants();
	unsigned long long seed = parser->getSeed();
	bool verbose = parser->hasVerbose();
	bool debug = parser->hasDebug();
	bool time = parser->hasTime();
	bool all_confs = parser->hasAllConfs();
	bool stability_only = parser->hasStabilityOnly();
	bool early_stop = parser->hasEarlyStop();
	bool states_only = parser->hasStatesOnly();

	// Timing variables
	std::chrono::time_point<std::chrono::high_resolution_clock> start_algo, stop_algo, start_overall;
	start_overall = std::chrono::high_resolution_clock::now();

	std::vector<double> rrc_vector_v = parser->getRRCVector();
	std::vector<int> state_change_matrix_v = parser->getStateChangeMatrix();
	std::vector<int> reactants_table_v = parser->getReactantsTable();
	std::vector<int> configuration_matrix_v = parser->getConfigurationMatrix();
	std::vector<double> propensity_matrix_v = parser->getPropensityMatrix();

	// Cast vector to raw data pointers
	double* rrc_vector = rrc_vector_v.data();
	int* state_change_matrix = state_change_matrix_v.data();
	int* configuration_matrix = configuration_matrix_v.data();
	int* reactants_table = reactants_table_v.data();
	double* propensity_matrix = propensity_matrix_v.data();

	int n = parser->getNumSpecies();
	int m = parser->getNumReactions();

	if (verbose) {
		std::cout << "Beginning simulation..." << std::endl << std::endl;
	}
	if (time) {
		start_algo = std::chrono::high_resolution_clock::now();
	}

	directMethod(state_change_matrix,
				 rrc_vector, 
				 configuration_matrix, 
				 propensity_matrix, 
				 reactants_table, 
				 sims, 
				 n, 
				 m, 
				 max_reactants, 
				 iters, 
				 verbose,
				 all_confs,
				 stability_only,
				 early_stop,
				 debug,
				 states_only,
				 rng,
				 seed,
				 partial_copy_indices);

	if (time) {
		stop_algo = std::chrono::high_resolution_clock::now();
	}

	if (time) {
		auto duration_algo = std::chrono::duration_cast<std::chrono::milliseconds>(stop_algo - start_algo);
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_algo - start_overall);

		std::cout << std::endl;
		std::cout << "~~~~~~~~~~~~~~~~~~~~Runtime Statistics~~~~~~~~~~~~~~~~~~~~" << std::endl;
		std::cout << "Simulation run time duration: " << duration_algo.count() / 1000.0 << " seconds." << std::endl;
		std::cout << "Total run time duration: " << duration.count() / 1000.0 << " seconds." << std::endl;
	}

	return 0;
}