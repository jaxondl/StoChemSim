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
    cout << ']\n';

    cout << "Start Props: [";
    for(int i = 0; i < my_parser.get_start_props().size(); i++)
        cout << my_parser.get_start_state().at(i) << ' ';
    cout << ']';



    return 0;
}