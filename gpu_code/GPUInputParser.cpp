//
// Created by Tarek on 10/16/2020.
//

#include "GPUInputParser.h"

// Reaction
GPUInputParser::Reaction::Reaction(double rrc, int *update_vector, int *reactant_coefs, string *rnames) {

}

double GPUInputParser::Reaction::calculate_propensity(unordered_map<string, int> species) {

}

// GPUInputParser

GPUInputParser::GPUInputParser() {

}

GPUInputParser::GPUInputParser(string fp) {

}

void GPUInputParser::add_reaction(double rrc, int *update_vector, int *reactant_coefs, string *names) {

}

void GPUInputParser::add_species(string name, int count) {

}

void GPUInputParser::process() {

}

int* GPUInputParser::get_start_state() {

}

int* GPUInputParser::get_start_props() {

}

GPUInputParser::Reaction* GPUInputParser::get_reactions() {

}

int** GPUInputParser::get_state_update_matrix() {

}