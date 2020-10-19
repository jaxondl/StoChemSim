//
// Created by Tarek on 10/16/2020.
//

#include "GPUInputParser.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
using namespace std;

// Reaction
GPUInputParser::Reaction::Reaction(double rrc, int *update_vector, int *rcoefs, string *rnames) {
    this->rrc = rrc;
    this->update_vector = update_vector;
    this->rcoefs = rcoefs;
    this->rnames = rnames;
}

double GPUInputParser::Reaction::calculate_propensity(unordered_map<string, int> species) {

}

// GPUInputParser

GPUInputParser::GPUInputParser() throws FileNotFoundException {
    ifstream infile("input/CRN.txt");
    if(infile.is_open()) {
        string line = "";

        // First line is # species
        getline(infile, line);
        stringstream stream(line);
        this.num_species = 0;
        stream >> num_species;

        // Second line is # reactions
        getline(infile, line);
        stringstream stream(line);
        this.num_reacts = 0;
        stream >> num_reacts;

        this.lines = new vector<string>(num_reacts, "");

        while(getline(infile, line)){
            vector.push_back(line)
        }
    }
    else {
        throw new FileNotFoundException("input/CRN.txt not found")
    }
}

GPUInputParser::GPUInputParser(string fp) {
    ifstream infile(fp);
    if(infile.is_open()) {
        string line = "";

        // First line is # species
        getline(infile, line);
        stringstream stream(line);
        this.num_species = 0;
        stream >> num_species;

        // Second line is # reactions
        getline(infile, line);
        stringstream stream(line);
        this.num_reacts = 0;
        stream >> num_reacts;

        this.lines = new vector<string>(2 + num_reacts, "");
        this.lines.push_back(to_string(num_species))
        this.lines.push_Back(to_string(num_reacts))

        while(getline(infile, line)){
            vector.push_back(line)
        }
    }
    else {
        throw new FileNotFoundException("input file not found")
    }
}

void GPUInputParser::add_reaction(double rrc, int *update_vector, int *reactant_coefs, string *names) {

}

void GPUInputParser::add_species(string name, int count) {

}

void GPUInputParser::process() {
    for(auto& line : this.lines) {
        for(auto& c : line) {
            if(isdigit(c)) {
                // decide if we're before or after the name, then collect total or number in reaction as appropriate
            }
            else if(c == "'") {
                // Get the name of the species character by character until we reach ' again
            }

        }
    }
}

vector<int> GPUInputParser::get_start_state() {

}

vector<int> GPUInputParser::get_start_props() {

}

GPUInputParser::Reaction* GPUInputParser::get_reactions() {

}

vector<vector<int>> GPUInputParser::get_state_update_matrix() {

}