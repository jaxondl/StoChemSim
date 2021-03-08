#include "directSSA/directMethodSSA.h"
#include "directSSA/reactionTree.h"
#include "common/decoder.h"
#include "common/dependencyGraph.h"
#include "common/inputVerifier.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream> 

using namespace std;

int main(int argc, char** argv) {
    inputVerifier *iv = new inputVerifier();
    string inputFilePath = argv[1]; // example: C:\\Users\\Isaac\\CLionProjects\\SeniorDesign\\CPU\\SampleInputs\\sample_input_SSA_file.txt
    string outputFilePath = argv[2]; // example: C:\\Users\\Isaac\\CLionProjects\\SeniorDesign\\CPU\\output.txt
    double endValue = stod(argv[3]);
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

        // Create SSA object and begin algorithm
        bool so = false;
        bool fo = false;
        bool ti = false;
        bool it = false; // iterations limit
        if(endValue <= 0)
            ti = true;
        else
            cout << "endValue is " << endValue << endl;
        for(int i = 4; i < argc; i++){
            string argument = argv[i];
            if (argument =="-so")
                so = true;
            else if (argument == "-fo")
                fo = true;
            else if (argument == "-it")
                it = true;
            
        }

        directMethodSSA *directSSA = new directMethodSSA(moleculeAmounts, reactionRates, reactantsVector, stateChangeVector, endValue, so, fo, ti, it);
        directSSA->start();

        // Print the results of the SSA - all states at all times, sequentially
        vector<double> allTimes = directSSA->getAllTimes();
        vector<vector<int>> allStates = directSSA->getAllStates();
        vector<int> currentState = directSSA->getCurrentState();
        double currentTime = directSSA->getCurrentTime();
        int currentIteration = directSSA->getCurrentIteration();

        ofstream outfile;
        outfile.open(outputFilePath);
        int iterationWidth = 15;
        char separator = ' ';
        if(fo){
            if(so) // final state only
                outfile << "Final State:" << endl;
            else { // final time and state only
                if(it){
                    outfile << "Final Iteration" << "\t\t\tFinal State" << endl;
                    outfile << left << setw(iterationWidth) << setfill(separator) << currentIteration << "\t\t\t";
                }
                else {
                    outfile << "Final Time" << "\t\t\tFinal State" << endl;
                    outfile << setprecision(5) << scientific << currentTime << "\t\t\t";
                }
            }
            for(int i = 0; i < currentState.size(); i++){
                    outfile << speciesList[i];
                    outfile << ": " << currentState[i];
                    outfile << "\t";
            }
        }
        else if(so){ // only states (NO TIME)
            outfile << "Iteration" << "\t\t\tState" << endl;
            for(int i = 0; i < allStates.size(); i++){
                outfile << left << setw(iterationWidth) << setfill(separator) << i << "\t\t";
                for(int j = 0; j < allStates[i].size(); j++){
                    outfile << speciesList[j];
                    outfile << ": " << allStates[i][j];
                    outfile << "\t";
                }
                outfile << endl;
            }
        }
        else{
            outfile << "Iteration" << "\t\t\tTime(s)" << "\t\t\tState" << endl;
            for(int i = 0; i < allTimes.size(); i++){
                outfile << left << setw(iterationWidth) << setfill(separator) << i << "\t\t" << setprecision(5) << scientific << allTimes[i] << "\t\t";
                for(int j = 0; j < allStates[i].size(); j++){
                    outfile << speciesList[j];
                    outfile << ": " << allStates[i][j];
                    outfile << "\t";
                }
                outfile << endl;
            }
        }
        outfile.close();
    }
    else {
        cout << "There were issue(s) with your input file." << endl;
    }
    return 0;
}
