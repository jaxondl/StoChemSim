#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

// Macro sets status to new status only if new status > status (ERROR > WARNING > ALERT)
#define priority(status, new_status) ((status==SAFE || (status == ALERT && new_status == ALERT) || (status != ERROR && new_status == WARNING)) ? new_status : status)

namespace po = boost::program_options;

enum arg_status_t { ERROR, WARNING, ALERT, SAFE };

// Create an options description; needs to be passed the addresses of variable storage locations
po::options_description create_opts_desc(int* sims, int* iters, unsigned long long* seed, std::string* infile, std::string* rng, std::vector<int>* partial_copy_indices);

// Create a positional options description
po::positional_options_description create_pos_opts_desc();

// Creat the variables map from the arguments and the descriptions
po::variables_map create_varmap(int argc, char* argv[], po::options_description desc, po::positional_options_description pos_desc);

void fill_implicit_vars(po::variables_map vm, bool* verbose, bool* time, bool* all_confs, bool* stability_only, bool* early_stop, bool* debug, bool* states_only);

// Validate the arguments stord in the variables map
arg_status_t validate_args(po::variables_map vm, int s, int n, int m);

// Validate input file seperately, so that parameters can be loaded
arg_status_t validate_infile(po::variables_map vm);

bool print_warning_message();

void print_error_message();

void print_safe_message();

#endif