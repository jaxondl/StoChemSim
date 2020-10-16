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
    unordered_map<string, int> species; // key is the species name, value is the tuple (index, count)
    Reaction* reactions; // All reactions in the CRN

    void add_reaction(double rrc, int* update_vector, int* reactant_coefs, string* names);
    void add_species(string name, int count);
public:
    GPUInputParser(); // assumes the file is called "CRN.txt" and is located in "../input/"
    GPUInputParser(string fp); // takes a string giving the location of the text file instead of assuming it's in input
    void process(); // Fills species and reactions vars by processing the input line by line
    int* get_species(); // returns a list of tuples of the structure (species name, count), ordered by index in species' values
    Reaction* get_reaction(); // returns a list of the reactions, ordered according to their occurence in the input file


};

#endif //CRN_SSA_WOLFRAM_PKG_GPUINPUTPARSER_H
