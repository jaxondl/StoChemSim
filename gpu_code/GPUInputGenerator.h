//
// Created by Tarek on 10/19/2020.
//
#include <string>
#include <vector>
using namespace std;

#ifndef CRN_SSA_WOLFRAM_PKG_GPUINPUTGENERATOR_H
#define CRN_SSA_WOLFRAM_PKG_GPUINPUTGENERATOR_H


class GPUInputGenerator {
private:
    int num_species;
    int num_reactions;
    int scaler;
    int max_coef;
    const vector<char> alpha = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V'};

    vector<string> generate_names(int n);
public:
    GPUInputGenerator(int num_species, int num_reactions, int scaler, int max_coef);
    void write_output(const string& fname);
};


#endif //CRN_SSA_WOLFRAM_PKG_GPUINPUTGENERATOR_H
