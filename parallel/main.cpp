//
// Created by Tarek on 10/19/2020.
//
#include "Parsing/GPUInputParser.h"
#include "Parsing/GPUInputGenerator.h"
#include <iostream>
using namespace std;

void test_gpu_input_parser(const string& fname) {
    GPUInputParser my_parser = GPUInputParser(fname);
    my_parser.process();

    cout << "Start State: [";
    for(int i = 0; i < my_parser.get_start_state().size(); i++)
        cout << my_parser.get_start_state().at(i) << ' ';
    cout << ']' << '\n';

    cout << "Start Props: [";
    for(int i = 0; i < my_parser.get_start_props().size(); i++)
        cout << my_parser.get_start_props().at(i) << ' ';
    cout << ']';

    cout << "\nState Update Matrix:\n";
    for(int i = 0; i < my_parser.get_num_reactions(); i++){
        cout << "[ ";
        for(int j = 0; j < my_parser.get_num_species(); j++){
            cout << my_parser.get_state_update_matrix()[i][j] << ' ';
        }
        cout << ']' << '\n';
    }

    cout << "\nIndex Keys:\n";
    unordered_map<string, int> myMap = my_parser.get_index_keys();
    for(const auto& elem : myMap) {
        std::cout << elem.first << " " << elem.second << " " << "\n";
    }
}

int test_gpu_input_generator(int num_species, int num_reactions, int scaler, int max_coef, const string& fname) {
    GPUInputGenerator my_generator(num_species, num_reactions, scaler, max_coef);
    my_generator.write_output(fname);
}

int main() {
    test_gpu_input_parser("customCRN.txt");
    //test_gpu_input_generator(100, 3000, 3, 5, "generatedCRN.txt");
    return 0;
}