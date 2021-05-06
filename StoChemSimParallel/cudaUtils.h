#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

// Write config array to binary file
void save_config(std::string fname, std::vector<int> config);

// Write times array to binary file
void save_times(std::string fname, std::vector<double> times);

// Check if stability flags are all true
bool is_stable(bool* stability_flags, int n);

#endif
