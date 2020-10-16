//
// Created by Tarek on 10/16/2020.
//
#include <unordered_map>
using namespace std;

#ifndef CRN_SSA_WOLFRAM_PKG_GPUINPUTPARSER_H
#define CRN_SSA_WOLFRAM_PKG_GPUINPUTPARSER_H

class GPUInputParser {
    class Reaction {
    private:
        double rrc; // rate reaction constant
        int* update_vector; // can be added to the count vector to update it if this reaction occurs
        int* rcoefs; // coefficients of the reactants for calculating propensity
        string* rnames; // Keys to get the counts of the reactants.
    public:
        Reaction(double rrc, int* update_vector, int* reactant_coefs, string* rnames);
        double calculate_propensity(unordered_map<string, int> species); // take the species map, reaction calculates its own propensity from it
    };

private:
    string* lines; // lines of the input text file
    unordered_map<string, int*> species; // key is the species name, value is the list [index, count]
    Reaction* reactions; // All reactions in the CRN

    void add_reaction(double rrc, int* update_vector, int* reactant_coefs, string* names);
    void add_species(string name, int count);
public:
    GPUInputParser(); // assumes the file is called "CRN.txt" and is located in "../input/"
    GPUInputParser(string fp); // takes a string giving the location of the text file instead of assuming it's in input
    void process(); // Fills species and reactions variables by processing the input line by line
    int* get_start_state(); // returns an array with the count of each chemical species ordered alphabetically (A-Z)
    int* get_start_props(); // returns an array with the starting propensity of each reaction, ordered by occurence in the input file.
    Reaction* get_reactions(); // returns an array containing all the reactions ordered by occurence in the input file.
    int** get_state_update_matrix(); // returns a matrix of the form [[-1, 1, -1, 0][2, -2, 0, 0]...], where each row updates the species count for a given reaction
};

#endif //CRN_SSA_WOLFRAM_PKG_GPUINPUTPARSER_H
