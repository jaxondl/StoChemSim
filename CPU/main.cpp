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
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\SampleInputs\\sample_deterministic_input_file_2.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_1.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_2.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_3.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_4.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\faulty_sample_input_SSA_file_5.txt");
    bool safeToRun = true;
    if (safeToRun) {
        decoder *inputDecoder = new decoder();
        inputDecoder->decode("C:\\Users\\Isaac\\CLionProjects\\SeniorDesign\\CPU\\SampleInputs\\sample_input_SSA_file.txt");
        vector<vector<pair<int, int>>> stateChangeVector = inputDecoder->getStateChangeVector();
        vector<vector<pair<int, int>>> reactantsVector = inputDecoder->getReactantVector();
        vector<double> reactionRates = inputDecoder->getkValueVector();
        vector<int> moleculeAmounts = inputDecoder->getPopulationSizes();
        algorithm *directSSA = new algorithm(moleculeAmounts, reactionRates, reactantsVector, stateChangeVector, 5.9);
        directSSA->start();
        //inputDecoder->decode("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\SampleInputs\\sample_deterministic_input_file_2.txt");
    }
//        /**
//        int numReactions = 2;
//        string moleculeTypes[6] = {"C", "B", "D", "E", "F", "A"};
//        int moleculeAmounts[6] = {314, 1235, 0, 0, 0, 2345};
//        vector <vector<pair<int, int>>> stateChangeArray(numReactions);
//        vector <vector<pair<int, int>>> reactantsArray(numReactions);
//        int kValues[2] = {3, 5};
//
//        pair<int, int> reactants1_1;
//        reactants1_1.first = 4;
//        reactants1_1.second = -1;
//        stateChangeArray[0].push_back(reactants1_1);
//        pair<int, int> reactants1_2;
//        reactants1_2.first = 4;
//        reactants1_2.second = -1;
//        stateChangeArray[0].push_back(reactants1_2);
//        pair<int, int> products1;
//        products1.first = 0;
//        products1.second = 1;
//        stateChangeArray[0].push_back((products1));
//
//        pair<int, int> reactants2_1;
//        reactants2_1.first = 3;
//        reactants2_1.second = -1;
//        stateChangeArray[1].push_back(reactants2_1);
//        pair<int, int> products2;
//        products2.first = 5;
//        products2.second = 1;
//        stateChangeArray[1].push_back(products2);
//
//        pair<int, int> rxn1_1;
//        rxn1_1.first = 4;
//        rxn1_1.second = 1;
//        reactantsArray[0].push_back(rxn1_1);
//        pair<int, int> rxn1_2;
//        rxn1_2.first = 1;
//        rxn1_2.second = 1;
//        reactantsArray[0].push_back(rxn1_2);
//
//        pair<int, int> rxn2_1;
//        rxn2_1.first = 3;
//        rxn2_1.second = 1;
//        reactantsArray[1].push_back(rxn2_1);
//
//        **/
//
//        cout << "Please input the method you would like to implement: " << endl;
//        int selection;
//        cin >> selection;
//
//        decoder decoderObject;
//        //dependency dependencyObject(moleculeTypes, moleculeAmounts, stateChangeArray, reactantsArray);
//
//
//        //tree test
//        string moleculeTypes[] = {"A","B","C","D","E"};
//        int moleculeAmounts[] = {5, 3, 3, 0, 0};
//        vector<pair<int, int>> ReactantsArray[3];
//        vector<pair<int,int>> reaction1;
//        pair<int,int> reactant1 = {1, 1};
//        pair<int,int> reactant2 = {2, 1};
//        reaction1.push_back(reactant1);
//        reaction1.push_back(reactant2);
//        ReactantsArray[0] = reaction1;
//
//        vector<pair<int,int>> reaction2;
//        pair<int,int> reactant3 = {1, 1};
//        pair<int,int> reactant4 = {2, 2};
//        reaction2.push_back(reactant3);
//        reaction2.push_back(reactant4);
//        ReactantsArray[1] = reaction2;
//
//        vector<pair<int,int>> reaction3;
//        pair<int,int> reactant5 = {1, 2};
//        pair<int,int> reactant6 = {2, 1};
//        reaction3.push_back(reactant5);
//        reaction3.push_back(reactant6);
//        ReactantsArray[2] = reaction3;
//
//        double rates[] = {1.25, 2, 2.35};
//
//        tree treeObj(moleculeAmounts, rates, ReactantsArray);
//        tree* tree_ptr = &treeObj;
//        //tree* tree = new tree(3, moleculeAmounts, rates, ReactantsArray);
//        moleculeAmounts[1] = 2;
//        tree_ptr->updatePropensity(1, rates[1], moleculeAmounts, ReactantsArray[1]);
//        tree_ptr->updatePropensity(2, rates[2], moleculeAmounts, ReactantsArray[2]);
//        for (int i = 0; i < 3; i++) {
//            cout << "At index " << i << endl;
//            cout << tree_ptr->ReactionTreeArray[i].propensity << endl;
//            cout << tree_ptr->ReactionTreeArray[i].leftChild << endl;
//            cout << tree_ptr->ReactionTreeArray[i].rightChild << endl;
//            cout << tree_ptr->ReactionTreeArray[i].leftSum << endl;
//            cout << tree_ptr->ReactionTreeArray[i].rightSum << endl;
//        }
//        cout << tree_ptr->searchForNode(0.48);
//
//
//
//
//        if (selection == 1) { //refer to the optimized SSA
//            algorithm algorithmObject;
//            //tree treeObject;
//
//            //algorithmObject.testFunction();
//            decoderObject.testFunction();
//            //dependencyObject.testFunction();
//            //treeObject.testFunction();
//        } else if (selection == 2) { //refer to tau-leaping, etc. (future reactions)
//            //dependencyObject.testFunction();
//            /** method 2 specific objects **/
//        }
//
//    } //this brace closes if(safeToRun)
    return 0;
}
