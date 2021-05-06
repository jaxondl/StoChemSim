#ifndef ARGPARSER_H
#define ARGPARSER_H

#pragma warning( disable : 4267)

#include <boost\program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "inputParser.h"

// Macro sets status to new status only if new status > status (ERROR > WARNING > ALERT)
#define priority(status, new_status) ((status==SAFE || (status == ALERT && new_status == ALERT) || (status != ERROR && new_status == WARNING)) ? new_status : status)
#define MAX_SIMS 65535 // This is the maximum dimension of the y/z dimension of a cuda block grid (2^31 - 1). TODO: Increase by using the z dimension for simulations as well.

namespace po = boost::program_options;

enum arg_status_t { ERROR, WARNING, ALERT, SAFE, HELP};

class ArgParser {
	private:
		po::options_description desc;
		po::positional_options_description pos_desc;
		po::variables_map vm;

		int argc;
		char** argv;

		std::vector<double> rrc_vector_v;
		std::vector<int> state_change_matrix_v;
		std::vector<int> reactants_table_v;
		std::vector<int> configuration_matrix_v;
		std::vector<double> propensity_matrix_v;
		int max_reactants;
		int m;
		int n;

		bool verbose = false;
		bool debug = false;
		bool time = false;
		bool all_confs = false;
		bool stability_only = false;
		bool early_stop = false;
		bool states_only = false;
		unsigned long long seed;
		int iters;
		int sims;
		std::vector<int> partial_copy_indices;
		std::string rng;
		std::string infile;

		// Create a positional options description
		po::positional_options_description create_pos_opts_desc();

		// Creat the variables map from the arguments and the descriptions
		po::variables_map create_varmap(int argc, char* argv[], po::options_description desc, po::positional_options_description pos_desc);

		void fill_implicit_vars(po::variables_map vm, bool* verbose, bool* time, bool* all_confs, bool* stability_only, bool* early_stop, bool* debug, bool* states_only);

		// Validate the arguments stord in the variables map
		arg_status_t validate_args(po::variables_map vm, int s, int n, int m);

		// Validate input file seperately, so that parameters can be loaded
		arg_status_t validate_infile(po::variables_map vm);

	public:

		// Constructs options description, positional options description, and variables map
		ArgParser(int argc, char* argv[]);

		// Reads and validates data from variables map, returning worst status code encountered
		arg_status_t parse();

		bool print_warning_message();

		void print_error_message();

		void print_safe_message();

		bool hasHelp();

		bool hasVerbose();

		bool hasDebug();

		bool hasTime();

		bool hasAllConfs();

		bool hasStabilityOnly();

		bool hasEarlyStop();

		bool hasStatesOnly();

		int getSims();

		int getIters();

		unsigned long long getSeed();

		std::string getInfile();

		std::string getRNG();

		std::vector<int> getPartialCopyIndices();
		
		std::vector<double> getRRCVector();

		std::vector<int> getStateChangeMatrix();

		std::vector<int> getReactantsTable();

		std::vector<int> getConfigurationMatrix();

		std::vector<double> getPropensityMatrix();

		int getMaxReactants();

		int getNumSpecies();

		int getNumReactions();
};

#endif