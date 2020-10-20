//
// Created by Tarek on 10/16/2020.
//

#include "GPUInputParser.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <utility>
using namespace std;

// Reaction
GPUInputParser::Reaction::Reaction(double rrc, vector<int> update_vector, vector<int> rcoefs, vector<int> rindices) {
    this->rrc = rrc;
    this->update_vector = move(update_vector);
    this->rcoefs = move(rcoefs);
    this->rindices = move(rindices);
}

double GPUInputParser::Reaction::calculate_propensity(vector<int> state) {
    double prop = rrc;
    for(int i = 0; i < rindices.size(); i++){
        int count = state.at(rindices.at(i));
        int coef = rcoefs.at(i);
        prop *= pow(count, coef);
    }
    return prop;
}

vector<int> GPUInputParser::Reaction::get_update_vector() {
    return update_vector;
}

// GPUInputParser

GPUInputParser::GPUInputParser() {
    ifstream infile;
    infile.open("../GPU/input/CRN.txt", ios_base::app|ios_base::in|ios_base::out);
    if(infile.is_open()) {
        string line;

        // First line is # species
        getline(infile, line);
        stringstream stream(line);
        num_species = 0;
        stream >> num_species;

        // Second line is # reactions
        getline(infile, line);
        stringstream stream2(line);
        num_reacts = 0;
        stream2 >> num_reacts;

        while(getline(infile, line)){
            lines.push_back(line);
        }
    }
    else {
        throw runtime_error("input/CRN.txt not found");
    }
}

GPUInputParser::GPUInputParser(const string& fname) {
    ifstream infile;
    infile.open("../GPU/input/" + fname, ios_base::app|ios_base::in|ios_base::out);
    if(infile.is_open()) {
        string line;
        num_species = 0;
        num_reacts = 0;

        // First line is # species
        getline(infile, line);
        stringstream stream(line);
        stream >> num_species;

        // Second line is # reactions
        getline(infile, line);
        stringstream stream2(line);
        stream2 >> num_reacts;

        while(getline(infile, line)){
            lines.push_back(line);
        }
    }
    else {
        throw runtime_error("input file " + fname + " not found");
    }
}

void GPUInputParser::process() {
    string top_del = ";"; // split into left, right, and rrc parts
    string mid_del = ","; // split into XX'ABC'YY components
    string low_del = "'"; // split into XX, ABC, YY tokens

    vector<int> tmp(num_species, 0);
    start_state = tmp;

    int index = 0;
    for(auto& line : lines) {
        double rrc = 0;
        vector<int> update_vector(num_species, 0);
        vector<int> rcoefs;
        vector<int> rindices;

        vector<string> reaction_parts = tokenize(line, top_del);
        string left = reaction_parts[0];
        string right = reaction_parts[1];

        stringstream stream(reaction_parts[2]);
        stream >> rrc;

        for(auto& component : tokenize(left, mid_del)){
            bool pre_name = true;
            int coef = 1;
            int count = -1;
            int rindex = -1;
            for(auto& token : tokenize(component, low_del)){
                // process the token as a species name
                if(!is_integer(token)){
                    pre_name = false;
                    // if the species hasn't been encountered, give it an index in the map then increment to get ready for next time.
                    if(index_keys.find(token) == index_keys.end()){
                        index_keys[token] = index;
                        index++;
                    }
                    rindex = index_keys[token];
                }
                // process the token as a coefficient
                else if(pre_name){
                    coef = 0;
                    stringstream stream2(token);
                    stream2 >> coef;
                }
                // process the token as a total count
                else{
                    count = 0;
                    stringstream stream3(token);
                    stream3 >> count;
                }
            }
            rcoefs.push_back(coef);
            rindices.push_back(rindex);
            if(count != -1){
                start_state.at(rindex) = count; // store count for first occurence
            }
            update_vector[rindex] -= coef; // subtract coef since it is reactant
        }
        for(auto& component : tokenize(right, mid_del)){
            bool pre_name = true;
            int coef = 1;
            int count = -1;
            int rindex = -1;
            for(auto& token : tokenize(component, low_del)){
                // process the token as a species name
                if(!is_integer(token)){
                    pre_name = false;
                    // if the species hasn't been encountered, give it an index in the map then increment to get ready for next time.
                    if(index_keys.find(token) == index_keys.end()){
                        index_keys[token] = index;
                        index++;
                    }
                    rindex = index_keys[token];
                }
                    // process the token as a coefficient
                else if(pre_name){
                    coef = 0;
                    stringstream stream2(token);
                    stream2 >> coef;
                }
                    // process the token as a total count
                else{
                    count = 0;
                    stringstream stream3(token);
                    stream3 >> count;
                }
            }
            if(count != -1){
                start_state[rindex] = count; // store count for first occurence
            }
            update_vector[rindex] += coef; // add coef since it is a product
        }

        GPUInputParser::Reaction react(rrc, update_vector, rcoefs, rindices);

        state_update_matrix.push_back(update_vector);
        start_props.push_back(react.calculate_propensity(start_state));
        reactions.push_back(react);
    }
}

unordered_map<string, int> GPUInputParser::get_index_keys(){
    return index_keys;
}

vector<int> GPUInputParser::get_start_state() {
    return start_state;
}

vector<double> GPUInputParser::get_start_props() {
    return start_props;
}

vector<GPUInputParser::Reaction> GPUInputParser::get_reactions() {
    return reactions;
}

vector<vector<int>> GPUInputParser::get_state_update_matrix() {
    return state_update_matrix;
}

vector<string> GPUInputParser::tokenize(string s, const string& delimiter) {
    size_t pos = 0;
    string token;
    vector<string> tokenized;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        tokenized.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    tokenized.push_back(s);
    return tokenized;
}

bool GPUInputParser::is_integer(string s){
    string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

int GPUInputParser::get_num_species() {
    return num_species;
}

int GPUInputParser::get_num_reactions() {
    return num_reacts;
}
