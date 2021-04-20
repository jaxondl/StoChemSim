#include "utils.h"

po::options_description create_opts_desc(int* sims, int* iters, unsigned long long* seed, std::string* infile, std::string* rng, std::vector<int>* partial_copy_indices) {
	// Options description for named flags with header text
	po::options_description desc("Available options given below. Please see the documentation for more extensive info.");

	desc.add_options()
		("help,h", "Provides more details on command line options.\n")
		("infile,f", po::value<std::string>(infile), "The location of the input file.\n")
		("sims,s", po::value<int>(sims), "The number of simulations to run.\n")
		("iters,i", po::value<int>(iters), "The number of iteartions to run all simulations for.\n")
		("verbose,v", "Log execution information to stdout.\n")
		("time,t", "If using the --verbose flag, record execution times as well.\n")
		("states-only", "Do not draw exponential random variables for time.\n")
		("all-confs,a", "Save all configurations; this may increase runtime noticeably.")
		("stability-only,o", "Only stop once stability has been reached in all sims. [WARNING]: This may cause the program not to terminate for certain CRNs.\n")
		("early-stop,e", "Stops after the number iterations specified or once stability is reached.\n")
		("rng", po::value<std::string>(rng)->default_value("MT19937"), "Select a random generator type, default is MT19937. Other options are: XORWOW, MRG32K3A, MTGP32, PHILOX-4X32-10, SOBOL32, scrambledSOBOL32, SOBOL64, scrambledSOBOL64.\n")
		("seed", po::value<unsigned long long>(seed)->default_value(42ULL), "Set the random seed to use with the random number generator. Must be an integer, default is 42.\n")
		("partial-copy,p", po::value<std::vector<int>>(partial_copy_indices)->multitoken(), "Only copy back a subset of the species for logging.\n")
		("debug,d", "Prints debug information to the console; use a small number (<=50) of simulations.\n")
		;

	return desc;
}

po::positional_options_description create_pos_opts_desc() {
	// Positional (unnamed) options, added in order
	po::positional_options_description pos_desc;
	pos_desc.add("infile", 1);
	pos_desc.add("sims", 1);
	pos_desc.add("iters", 1);

	return pos_desc;
}

po::variables_map create_varmap(int argc, char* argv[], po::options_description desc, po::positional_options_description pos_desc) {
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
	po::notify(vm);

	return vm;
}

void fill_implicit_vars(po::variables_map vm, bool* verbose, bool* time, bool* all_confs, bool* stability_only, bool* early_stop, bool* debug, bool* states_only) {
	if (vm.count("verbose")) {
		*verbose = true;
	}
	if (vm.count("time")) {
		*time = true;
	}
	if (vm.count("all-confs")) {
		*all_confs = true;
	}
	if (vm.count("stability-only")) {
		*stability_only = true;
	}
	if (vm.count("early-stop")) {
		*early_stop = true;
	}
	if (vm.count("debug")) {
		*debug = true;
	}
	if (vm.count("states-only")) {
		*states_only = true;
	}
}

arg_status_t validate_args(po::variables_map vm, int s, int n, int m) {
	arg_status_t status = SAFE;

	if (vm.count("sims")) {
		if (s != -1) {
			std::cout << "[ALERT]: Number of simulations is specified in the input file and by command; the command will be used instead." << std::endl;
			s = vm["sims"].as<int>();
			status = priority(status, ALERT);
		}
	}
	else if (s == -1) {
		std::cout << "[ERROR]: Number of simulations not specified; must be in either input file or given by command line." << std::endl;
		status = ERROR;
	}

	if (vm.count("iters")) {
		if (vm.count("stability-only")) {
			std::cout << "[ALERT]: Flag --stability-only set, the number of iterations given will not be used." << std::endl;
			status = priority(status, ALERT);
		}
	}
	else if (!vm.count("stability-only")) {
		std::cout << "[ERROR]: Flag --stability-only must be set if the number of iterations is not specified." << std::endl;
		status = ERROR;
	}

	if (vm.count("time") && !vm.count("verbose")) {
		std::cout << "[WARNING]: Timing information will not be displayed unless --verbose is also set." << std::endl;
		status = priority(status, WARNING);
	}

	if (vm.count("all-confs")) {
		double conf_size;
		if (vm.count("partial-copy")) {
			int species = vm["partial-copy"].as<std::vector<int>>().size();
			conf_size = s * species * sizeof(int) / 1000000000.0; // number of GBs in a partial configuration matrix
		}
		else {
			conf_size = s * n * sizeof(int) / 1000000000.0; // number of GBs in a full configuration matrix
		}

		if (vm.count("stability-only")) {
			std::cout << "[ERROR]: --stability-only and all-confs both set; this is not allowed since it may generate unbounded amounts of data." << std::endl;
			status = ERROR;
		}
		else if (vm.count("iters") && vm["iters"].as<int>() * conf_size > 1) {
			double size = vm["iters"].as<int>() * conf_size;
			std::cout.precision(3);
			std::cout << "[WARNING]: Input may cause up to ~" << std::fixed << size << " gigabytes of data to be generated; ensure you have adequate disc space available." << std::endl;
			status = WARNING;
		}
	}

	if (vm.count("stability-only")) {
		std::cout << "[ALERT]: The flag --stability-only may cause the program to run forever for some CRNs." << std::endl;
		status = ALERT;
	}

	if (vm.count("rng")) {
		std::string rng = vm["rng"].as<std::string>();
		if (rng != "XORWOW" &&
			rng != "MRG32K3A" &&
			rng != "MTGP32" &&
			rng != "PHILOX-4X32-10" &&
			rng != "SOBOL32" &&
			rng != "SOBOL64" &&
			rng != "scrambledSOBOL32" &&
			rng != "scrambledSOBOL64" &&
			rng != "MT19937") 
		{
			std::cout << "[ERROR]: Random generator type not recognized; see help for a list of options." << std::endl;
			status = ERROR;
		}
	}

	if (vm.count("partial-copy")) {
		if (!vm.count("all-confs")) {
			std::cout << "[WARNING]: Partial copy will have no effect if --all-confs is not set." << std::endl;
			status = WARNING;
		}
	}

	return status;
}

arg_status_t validate_infile(po::variables_map vm) {
	arg_status_t status;

	if (vm.count("infile")) {
		std::string fpath = vm["infile"].as<std::string>();
		std::ifstream infile;
		infile.open(fpath);
		if (!infile) {
			std::cout << "[ERROR]: Input file does not exist." << std::endl;
			status = ERROR;
		}
	}
	else {
		std::cout << "[ERROR]: Input file not specified." << std::endl;
		status = ERROR;
	}

	return status;
}

bool print_warning_message() {
	std::string choice;

	std::cout << std::endl;
	std::cout << "[WARNING]: There are warnings present. Do you still want to continue? (y/n): ";
	std::cin >> choice;
	std::cout << std::endl;

	if (choice != "y") {
		std::cout << "Operation aborted." << std::endl;
		return false;
	}
	return true;
}

void print_error_message() {
	std::cout << std::endl;
	std::cout << "[ERROR]: Invalid arguments must be fixed before execution." << std::endl;
}

void print_safe_message() {
	std::cout << std::endl;
	std::cout << "[SAFE] No warnings or errors encountered during command verification." << std::endl;
}