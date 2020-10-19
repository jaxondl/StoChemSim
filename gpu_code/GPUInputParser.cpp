//
// Created by Tarek on 10/16/2020.
//

#include "GPUInputParser.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <cmath>
using namespace std;

// Reaction
GPUInputParser::Reaction::Reaction(double rrc, int *update_vector, int *rcoefs, string *rindices) {
    this->rrc = rrc;
    this->update_vector = update_vector;
    this->rcoefs = rcoefs;
    this->rindices = rindices;
}

double GPUInputParser::Reaction::calculate_propensity(vector<int> start_state) {
    double prop = rrc;
    for(i = 0; i < this.rindices.size(); i++){
        count = start_state[rindices[i]];
        coef = rccoefs[i];
        prop *= pow(count, coef)
    }
    return prop;
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
        this.num_species = 0;
        this.num_reacts = 0;

        // First line is # species
        getline(infile, line);
        stringstream stream(line);
        stream >> num_species;

        // Second line is # reactions
        getline(infile, line);
        stringstream stream(line);
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

void GPUInputParser::process() {
    string top_del = ";"; // split into left, right, and rrc parts
    string mid_del = ","; // split into XX'ABC'YY components
    string low_del = "'"; // split into XX, ABC, YY tokens

    int index = 0;
    for(auto& line : this.lines) {
        double rrc = 0;
        vector<int> update_vector(this.num_species, 0);
        vector<int> rcoefs;
        vector<int> rindices;

        vector<string> reaction_parts = this.tokenize(line, top_del);
        left = reaction_parts[0];
        right = reaction_parts[1];

        stringstream stream(reaction_parts[2]);
        stream >> rrc;

        for(auto& component : this.tokenize(left, mid_del)){
            bool pre_name = true;
            int coef = 1;
            int count = NAN;
            int rindex = NAN;
            for(auto& token : this.tokenize(component, low_del)){
                // process the token as a species name
                if(!is_integer(token)){
                    pre_name = false;
                    // if the species hasn't been encountered, give it an index in the map then increment to get ready for next time.
                    if(this.index_keys.find(token) == this.index_keys.end()){
                        this.index_keys[token] = index;
                        index++;
                    }
                    rindex = this.index_keys[token];
                }
                // process the token as a coefficient
                else if(pre_name){
                    coef = 0;
                    stringstream stream(token);
                    stream >> coef;
                }
                // process the token as a total count
                else{
                    stringstream stream(token);
                    count = 0;
                    stream >> count;
                }
            }
            rcoefs.push_back(coef);
            rindices.push_back(rindex);
            if(!isnan(count)){
                this.species[rindex] = count; // store count for first occurence
            }
            update_vector[rindex] -= coef; // subtract coef since it is reactant
        }
        for(auto& component : this.tokenize(right, mid_del)){
            bool pre_name = true;
            int coef = 1;
            int count = NAN;
            int rindex = NAN;
            for(auto& token : this.tokenize(component, low_del)){
                // process the token as a species name
                if(!is_integer(token)){
                    pre_name = false;
                    // if the species hasn't been encountered, give it an index in the map then increment to get ready for next time.
                    if(this.index_keys.find(token) == this.index_keys.end()){
                        this.index_keys[token] = index;
                        index++;
                    }
                    rindex = this.index_keys[token];
                }
                    // process the token as a coefficient
                else if(pre_name){
                    coef = 0;
                    stringstream stream(token);
                    stream >> coef;
                }
                    // process the token as a total count
                else{
                    stringstream stream(token);
                    count = 0;
                    stream >> count;
                }
            }
            if(!isnan(count)){
                this.species[rindex] = count; // store count for first occurence
            }
            update_vector[rindex] += coef; // add coef since it is a product
        }

        GPUInputParser::Reaction react(rrc, update_vector, rcoefs, rindices);

        state_update_matrix.push_back(update_vector)
        start_props.push_back(react.calculate_propensity())
        reactions.push_back(react);
    }
}

unordered_map<string, int> GPUInputParser::get_index_keys(){
    return index_keys;
}

vector<int> GPUInputParser::get_start_state() {
    return start_state;
}

vector<int> GPUInputParser::get_start_props() {
    return start_props;
}

vector<GPUInputParser::Reaction> GPUInputParser::get_reactions() {
    return reactions;
}

vector<vector<int>> GPUInputParser::get_state_update_matrix() {
    vector<vector<int>> state_update_matrix;
}

vector<string> GPUInputParser::tokenize(string s, string delimiter) {
    size_t pos = 0;
    string token;
    vector<string> tokenized;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        tokenized.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    return tokenized
}

bool GPUInputParser::is_integer(string s){
    string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}