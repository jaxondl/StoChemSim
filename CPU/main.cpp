#include "algorithm.h"
#include "decoder.h"
#include "dependency.h"
#include "tree.h"
#include "inputVerifier.h"

#include <iostream>
#include <string>

using namespace std;

int main() {

    //inputVerifier *iv = new inputVerifier();
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\sample_input_SSA_file.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\sample_deterministic_input_file_1.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_1.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_2.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_3.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_4.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_5.txt");

    //if (safeToRun) {
        int numReactions = 2;
        string moleculeTypes[6] = {"C", "B", "D", "E", "F", "A"};
        int moleculeAmounts[6] = {314, 1235, 0, 0, 0, 2345};
        vector <vector<pair<int, int>>> stateChangeArray(numReactions);
        vector <vector<pair<int, int>>> reactantsArray(numReactions);
        int kValues[2] = {3, 5};

        pair<int, int> reactants1_1;
        reactants1_1.first = 4;
        reactants1_1.second = -1;
        stateChangeArray[0].push_back(reactants1_1);
        pair<int, int> reactants1_2;
        reactants1_2.first = 4;
        reactants1_2.second = -1;
        stateChangeArray[0].push_back(reactants1_2);
        pair<int, int> products1;
        products1.first = 0;
        products1.second = 1;
        stateChangeArray[0].push_back((products1));

        pair<int, int> reactants2_1;
        reactants2_1.first = 3;
        reactants2_1.second = -1;
        stateChangeArray[1].push_back(reactants2_1);
        pair<int, int> products2;
        products2.first = 5;
        products2.second = 1;
        stateChangeArray[1].push_back(products2);

        pair<int, int> rxn1_1;
        rxn1_1.first = 4;
        rxn1_1.second = 1;
        reactantsArray[0].push_back(rxn1_1);
        pair<int, int> rxn1_2;
        rxn1_2.first = 1;
        rxn1_2.second = 1;
        reactantsArray[0].push_back(rxn1_2);

        pair<int, int> rxn2_1;
        rxn2_1.first = 3;
        rxn2_1.second = 1;
        reactantsArray[1].push_back(rxn2_1);

        cout << "Please input the method you would like to implement: " << endl;
        int selection;
        cin >> selection;

        decoder decoderObject;
        dependency dependencyObject(moleculeTypes, moleculeAmounts, stateChangeArray, reactantsArray);

        if (selection == 1) { //refer to the optimized SSA
            algorithm algorithmObject;
            tree treeObject;

            algorithmObject.testFunction();
            decoderObject.testFunction();
            dependencyObject.testFunction();
            treeObject.testFunction();
        } else if (selection == 2) { //refer to tau-leaping, etc. (future reactions)
            dependencyObject.testFunction();
            /** method 2 specific objects **/
        }
    //}
    return 0;
}
