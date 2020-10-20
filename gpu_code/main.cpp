//
// Created by Tarek on 10/19/2020.
//
#include "GPUInputParser.h"
#include <iostream>
using namespace std;


int main() {
    GPUInputParser my_parser = GPUInputParser();
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

    return 0;
}