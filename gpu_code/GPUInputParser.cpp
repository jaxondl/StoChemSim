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

void GPUInputParser::add_reaction(double rrc, int *update_vector, int *reactant_coefs, string *names) {

}

void GPUInputParser::add_species(string name, int count) {

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
        vector<int> rnames;

        vector<string> reaction_parts = this.tokenize(line, top_del);
        left = reaction_parts[0];
        right = reaction_parts[1];

        stringstream stream(reaction_parts[2]);
        stream >> rrc;

        for(auto& component : this.tokenize(left, mid_del)){
            bool pre_name = true;
            int coef = 1;
            int count = NAN;
            string name = "";
            for(auto& token : this.tokenize(component, low_del)){
                if(!is_integer(token)){
                    pre_name = false;
                    // if the species hasn't been encountered, give it an index then increment to get ready for next time.
                    if(this.index_keys.find(token) == this.index_keys.end()){
                        this.index_keys[token] = index;
                        index++;
                    }
                    name = token;
                }
                else if(pre_name){
                    // process it as a coefficient
                    coef = 0;
                    stringstream stream(token);
                    stream >> coef;
                }
                else{
                    // process it as a total count
                    stringstream stream(token);
                    count = 0;
                    stream >> count;
                }
            }
            rcoefs.push_back(coef);
            rnames.push_back(name);
            if(!isnan(count)){
                this.
            }
        }
        for(auto& component : this.tokenize(right, mid_del)){
            for(auto& token : this.tokenize(component, low_del)){

            }
        }
    }
}

unordered_map<string, int> GPUInputParser::get_index_keys(){
    
}

vector<int> GPUInputParser::get_start_state() {

}

vector<int> GPUInputParser::get_start_props() {

}

GPUInputParser::Reaction* GPUInputParser::get_reactions() {

}

vector<vector<int>> GPUInputParser::get_state_update_matrix() {

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