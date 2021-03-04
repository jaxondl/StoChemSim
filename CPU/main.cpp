#include "directMethodSSA.h"
#include "decoder.h"
#include "dependencyGraph.h"
#include "reactionTree.h"
#include "inputVerifier.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream> 

using namespace std;

int main(int argc, char** argv) {
    inputVerifier *iv = new inputVerifier();
    string inputFilePath = argv[1]; // example: C:\\Users\\Isaac\\CLionProjects\\SeniorDesign\\crn-ssa-wolfram-pkg\\CPU\\SampleInputs\\sample_input_SSA_file.txt
    string outputFilePath = argv[2]; // example: C:\\Users\\Isaac\\CLionProjects\\SeniorDesign\\crn-ssa-wolfram-pkg\\CPU\\output.txt
    bool safeToRun = iv->verifyFile(inputFilePath);
    if (safeToRun) {
        // Create decoder object and get all needed data structures
        decoder *inputDecoder = new decoder();
        inputDecoder->decode(inputFilePath);
        vector<vector<pair<int, int>>> stateChangeVector = inputDecoder->getStateChangeVector();
        vector<vector<pair<int, int>>> reactantsVector = inputDecoder->getReactantVector();
        vector<double> reactionRates = inputDecoder->getkValueVector();
        vector<int> moleculeAmounts = inputDecoder->getPopulationSizes();
        vector<string> speciesList = inputDecoder->getListOfSpecies();
        double tEnd = inputDecoder->getTEnd();

        // Create SSA object and begin algorithm
        directMethodSSA *directSSA = new directMethodSSA(moleculeAmounts, reactionRates, reactantsVector, stateChangeVector, tEnd);
        directSSA->start();

        // Print the results of the SSA - all states at all times, sequentially
        vector<double> allTimes = directSSA->getAllTimes();
        vector<vector<int>> allStates = directSSA->getAllStates();

        ofstream outfile;
        outfile.open(outputFilePath);
        outfile << "Iteration" << "\t\tTime(s) " << "\t\t\tState" << endl;
        for(int i = 0; i < allTimes.size(); i++){
            outfile << i << "\t\t\t\t" << setprecision(5) << scientific << allTimes[i] << "\t\t\t";
            for(int j = 0; j < allStates[i].size(); j++){
                outfile << speciesList[j];
                outfile << ": " << allStates[i][j];
                outfile << "\t";
            }
            outfile << endl;
        }
        outfile.close();
    }
    else {
        cout << "There were issue with your input file." << endl;
    }
    return 0;
}