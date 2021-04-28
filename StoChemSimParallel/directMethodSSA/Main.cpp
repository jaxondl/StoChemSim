#include "inputParser.h"
#include "DirectMethod.cuh"
#include "utils.h"
#include "cudaUtils.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost\program_options.hpp>

#define MAX_SIMS 65535 // This is the maximum dimension of the y/z dimension of a cuda block grid (2^31 - 1). TODO: Increase by using the z dimension for simulations as well.

namespace po = boost::program_options;

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
	// Command option variables
	std::string infile;
	std::string rng;
	std::vector<int> partial_copy_indices;
	int sims;
	int iters;
	unsigned long long seed;
	bool verbose = false;
	bool debug = false;
	bool time = false;
	bool all_confs = false;
	bool stability_only = false;
	bool early_stop = false;
	bool states_only = false;

	// Timing variables
	std::chrono::time_point<std::chrono::high_resolution_clock> start_algo, stop_algo, start_overall;
	start_overall = std::chrono::high_resolution_clock::now();

	// Perform command line parsing
	po::options_description desc = create_opts_desc(&sims, &iters, &seed, &infile, &rng, &partial_copy_indices);
	po::positional_options_description pos_desc = create_pos_opts_desc();
	po::variables_map vm = create_varmap(argc, argv, desc, pos_desc);

	// If help flag is set, output description and quit.
	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 0;
	}

	fill_implicit_vars(vm, &verbose, &time, &all_confs, &stability_only, &early_stop, &debug, &states_only);

	// Output banner
	if (verbose) {
		std::ifstream banner_file;
		banner_file.open("banner_logo.txt");

		std::string line;
		while (std::getline(banner_file, line)) {
			std::cout << line << std::endl;
		}
		std::cout << std::endl;

		std::cout << "Decoding input file... " << std::endl;
	}

	arg_status_t infile_status = validate_infile(vm);
	if (infile_status == ERROR) {
		print_error_message();
		return 0;
	}

	gpuDecoderPrototype* decoder = new gpuDecoderPrototype();
	bool s_found = decoder->decode(infile);

	int s = decoder->getNumSimulations();
	int m = decoder->getNumReactions();
	int n = decoder->getNumSpecies();

	if (verbose) {
		std::cout << "Validating commands... " << std::endl << std::endl;
	}

	// Will generate a log of errors and warnings. Execution halts in the presence of errors, 
	// but prompts user to continue if there are warnings
	arg_status_t status = validate_args(vm, s, n, m);

	// if sims was input, use it instead of the input file
	if (sims) {
		s = sims;
	}

	if (s > MAX_SIMS) {
		std::cout << "[ERROR]: Maximum number of simulations exceeded. Only 2^31 - 1 = 65535 simulations are allowed." << std::endl;
		status = ERROR;
	}

	if (status == ERROR) {
		print_error_message();
		return 0;
	}
	else if (status == WARNING) {
		bool proceed = print_warning_message();
		if (!proceed) { return 0; }
	}
	else if (status == ALERT || status == SAFE) {
		print_safe_message();
	}
	std::cout << std::endl;

	// time should only be true if verbose is true
	if (!verbose) {
		time = false;
	}

	// If we are not saving the configs, do not transfer the partials back since there is no point
	if (!all_confs) {
		partial_copy_indices.clear();
	}

	std::vector<double> rrc_vector_v = decoder->getRRCVector();
	std::vector<int> state_change_matrix_v = decoder->getStateChangeMatrix();
	std::vector<int> reactants_table_v = decoder->getReactantsTableVector();
	std::vector<int> configuration_matrix_v;
	std::vector<double> propensity_matrix_v;

	// The decoder uses the number of simulations to copy the initial configuration and propensity values s times, 
	// so if the command line is used instead then this must be passed to the decoder.
	configuration_matrix_v = decoder->getConfigurationMatrix(s);
	propensity_matrix_v = decoder->getPropensityMatrix(s);

	int max_reactants = reactants_table_v.size() / (2 * m); // reactants table is of shape max_reacants * 2 * m

	// Cast vector to raw data pointers
	double* rrc_vector = rrc_vector_v.data();
	int* state_change_matrix = state_change_matrix_v.data();
	int* configuration_matrix = configuration_matrix_v.data();
	int* reactants_table = reactants_table_v.data();
	double* propensity_matrix = propensity_matrix_v.data();

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
				 s, 
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