//
// Created by Tarek on 10/19/2020.
//

#include "GPUInputGenerator.h"
#include <fstream>
#include <unordered_map>
#include <cmath>
using namespace std;

GPUInputGenerator::GPUInputGenerator(int num_species, int num_reactions, int scaler, int max_coef) {
    this->num_species = num_species;
    this->num_reactions = num_reactions;
    this->scaler = scaler;
    this->max_coef = max_coef;
}

void GPUInputGenerator::write_output(const string& fname) {
    ofstream outfile;
    outfile.open("../gpu_code/output/" + fname, ios_base::app|ios_base::in|ios_base::out);
    if(outfile.is_open()){
        outfile << num_species << '\n';
        outfile << num_reactions << '\n';
        vector<string> names = generate_names(num_species);
        
        unordered_map<string, bool> encounter_map;
        bool all_names_used = false;
        int nindex = 0;
        for(int i=0; i < num_reactions; i++){
            int num_reactants = rand() % 5 + 1;
            int num_products = rand() % 3 + 1;
            double rrc = (rand() % 10 + 1.0) / 10.0;

            unordered_map<string, bool> species_map;
            for(int j = 0; j < num_reactants; j++){
                int coef = rand() % max_coef + 1;

                // If all names have been used, then use whichever. Otherwise generate them in order.
                // TODO: Better system for ensuring all names get used instead of just in-order output.
                string name;
                if(all_names_used){
                    int name_index = rand() % names.size();
                    name = names.at(name_index);
                }
                else{
                    name = names.at(nindex);
                    nindex++;
                    if(nindex == names.size()){
                        all_names_used = true;
                    }
                }

                // If the name has been used in the reaction, go back an iteration and re-generate the name.
                if(species_map.find(name) != species_map.end()){
                    j--;
                    continue;
                }
                species_map[name] = true;

                // If we have not encountered the species, we give it a count.
                int count = -1;
                if(encounter_map.find(name) == encounter_map.end()){
                    count = pow((rand() % 10 + 1), scaler);
                    encounter_map[name] = true;
                }

                if(coef != 1){
                    outfile << coef << "'";
                }

                outfile << name;

                if(count != -1){
                    outfile << "'" << count;
                }

                if(j != num_reactants - 1){
                    outfile << ",";
                }
            }
            outfile << ";";
            species_map.clear();
            for(int j = 0; j < num_products; j++){
                int coef = rand() % max_coef + 1;

                string name;
                if(all_names_used){
                    int name_index = rand() % names.size();
                    name = names.at(name_index);
                }
                else{
                    name = names.at(nindex);
                    nindex++;
                    if(nindex == names.size()){
                        all_names_used = true;
                    }
                }

                // If the name has been used in the reaction, go back an iteration and re-generate the name.
                if(species_map.find(name) != species_map.end()){
                    j--;
                    continue;
                }
                species_map[name] = true;

                // If we have not encountered the species, we give it a count.
                int count = -1;
                if(encounter_map.find(name) == encounter_map.end()){
                    count = pow((rand() % 10 + 1), scaler);
                    encounter_map[name] = true;
                }

                if(coef != 1){
                    outfile << coef << "'";
                }

                outfile << name;

                if(count != -1){
                    outfile << "'" << count;
                }

                if(j != num_products - 1){
                    outfile << ",";
                }
            }
            outfile << ";" << rrc;

            if(i != num_reactions  - 1){
                outfile << '\n';
            }
        }
    }
    else{
        throw runtime_error("couldn't open outfile " + fname);
    }

    outfile.close();
}

vector<string> GPUInputGenerator::generate_names(int n) {
    vector<string> names;
    int letters = 1;
    for(int i = 0; i < n; i++){
        string name;
        for(int l = 0; l < letters; l++){
            name += alpha[i % 26];
        }
        if(alpha[i % 26] == "Z"){
            letters++;
        }
        names.push_back(name);
    }
    return names;
}
