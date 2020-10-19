//
// Created by Tarek on 10/16/2020.
//
#include <unordered_map>
#include <string>
#include <vector>
using namespace std;

#ifndef CRN_SSA_WOLFRAM_PKG_GPUINPUTPARSER_H
#define CRN_SSA_WOLFRAM_PKG_GPUINPUTPARSER_H

class GPUInputParser {
    class Reaction {
    private:
        double rrc; // rate reaction constant
        vector<int> update_vector; // can be added to the count vector to update it if this reaction occurs
        vector<int> rcoefs; // coefficients of the reactants for calculating propensity
        vector<string> rindices; // Keys to get the counts of the reactants.
    public:
        Reaction(double rrc, vector<int> update_vector, vector<int> rcoefs, vector<string> rindices);
        double calculate_propensity(vector<int> start_state); // take the species map, reaction calculates its own propensity from it
    };

private:
    vector<string> lines; // excludes the first two lines (they are stored in num_species & num_reacts respectively)
    vector<int> start_state; // vector showing the count of each species indexed as given by index_keys
    vector<double> start_props;
    unordered_map<string, int> index_keys; // key is the species name, value is its numeric index
    vector<Reaction> reactions; // All reactions in the CRN
    vector<vector<int>> state_update_matrix;
    int num_species;
    int num_reacts;

    static vector<string> tokenize(string s, string delimiter);
    static bool is_integer(string s);
public:
    GPUInputParser(); // assumes the file is called "CRN.txt" and is located in "../input/"
    GPUInputParser(string fp); // takes a string giving the location of the text file instead of assuming it's in input
    void process(); // Fills species and reactions variables by processing the input line by line
    unordered_map<string, int> get_index_keys();
    vector<int> get_start_state(); // returns an array with the count of each chemical species ordered as they occur in the input
    vector<int> get_start_props(); // returns an array with the starting propensity of each reaction, ordered by occurence in the input file.
    vector<Reaction> get_reactions(); // returns an array containing all the reactions ordered by occurence in the input file.
    vector<vector<int>> get_state_update_matrix(); // returns a matrix of the form [[-1, 1, -1, 0][2, -2, 0, 0]...], where each row updates the species count for a given reaction
};

#endif //CRN_SSA_WOLFRAM_PKG_GPUINPUTPARSER_H
