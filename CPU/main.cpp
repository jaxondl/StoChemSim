#include "directMethodSSA.h"
#include "decoder.h"
#include "dependencyGraph.h"
#include "reactionTree.h"
#include "inputVerifier.h"

#include <iostream>
#include <string>

using namespace std;

int main() {
    inputVerifier *iv = new inputVerifier();
    bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\SampleInputs\\sample_input_SSA_file.txt");
    //bool safeToRun = true;
    if (safeToRun) {
        // Create decoder object and get all needed data structures
        decoder *inputDecoder = new decoder();
        //inputDecoder->decode("/Users/vidursinha/Desktop/Senior-Design/crn-ssa-wolfram-pkg/CPU/SampleInputs/sample_input_SSA_file.txt");
        inputDecoder->decode("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\CPU\\SampleInputs\\sample_input_SSA_file.txt");
        vector<vector<pair<int, int>>> stateChangeVector = inputDecoder->getStateChangeVector();
        vector<vector<pair<int, int>>> reactantsVector = inputDecoder->getReactantVector();
        vector<double> reactionRates = inputDecoder->getkValueVector();
        vector<int> moleculeAmounts = inputDecoder->getPopulationSizes();
        double tEnd = inputDecoder->getTEnd();

        // Create SSA object and begin algorithm
        directMethodSSA *directSSA = new directMethodSSA(moleculeAmounts, reactionRates, reactantsVector, stateChangeVector, tEnd);
        directSSA->start();

        // Print the results of the SSA - all states at all times, sequentially
        vector<double> allTimes = directSSA->getAllTimes();
        vector<vector<int>> allStates = directSSA->getAllStates();
        for(int i = 0; i < allTimes.size(); i++){
            cout << "Time: " << allTimes[i] << " State:";
            for(int stateNum: allStates[i]){
                cout << " " << stateNum; 
            }
            cout << endl;
         }
    }
    else {
        cout << "There were issue with your input file" << endl;
    }
    return 0;
}