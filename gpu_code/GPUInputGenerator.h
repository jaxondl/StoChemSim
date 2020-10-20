//
// Created by Tarek on 10/19/2020.
//
#include <string>
#include <vector>
#include <windows.h>
using namespace std;

#ifndef CRN_SSA_WOLFRAM_PKG_GPUINPUTGENERATOR_H
#define CRN_SSA_WOLFRAM_PKG_GPUINPUTGENERATOR_H


class GPUInputGenerator {
private:
    int num_species;
    int num_reactions;
    int scaler;
    int max_coef;
    const vector<string> alpha = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};

    vector<string> generate_names(int n);
public:
    GPUInputGenerator(int num_species, int num_reactions, int scaler, int max_coef); // num_reactions >= num_species or else num_species be used properly.
    void write_output(const string& fname); // Writes the output file according to the specifications given in the constructor.
};


#endif //CRN_SSA_WOLFRAM_PKG_GPUINPUTGENERATOR_H
