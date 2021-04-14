#include "common/decoder.h"
#include "common/inputVerifier.h"
#include "boundedTauLeaping/boundedTauLeaping.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {
    inputVerifier *iv = new inputVerifier(); // create input verifier for validating input file
    string inputFilePath = argv[1]; // example: C:\\Users\\Isaac\\CLionProjects\\SeniorDesign\\crn-ssa-wolfram-pkg\\sequential\\inputs\\sample_input_SSA_file.txt
    string outputFilePath = argv[2]; // example: C:\\Users\\Isaac\\CLionProjects\\SeniorDesign\\crn-ssa-wolfram-pkg\\sequential\\output.txt
    double endValue = stod(argv[3]); // a double used to specify the end time or iteration, depending on the user's flags
    bool safeToRun = iv->verifyFile(inputFilePath); // the input verifier will indicate whether the input file is valid and safe to run
    if (safeToRun) {
        // Create decoder object and get all needed data structures
        decoder *inputDecoder = new decoder();
        inputDecoder->decode(inputFilePath);
        vector<vector<pair<int, int> > > stateChangeVector = inputDecoder->getStateChangeVector();
        vector<vector<pair<int, int> > > reactantsVector = inputDecoder->getReactantVector();
        vector<double> reactionRates = inputDecoder->getkValueVector();
        vector<int> moleculeAmounts = inputDecoder->getPopulationSizes();
        vector<string> speciesList = inputDecoder->getListOfSpecies();

        // begin setting flags from the command line arguments, with every flag by default set to false
        //bool so = false; // states only flag
        bool fo = false; // final only flag
        bool ti = false; // infinite time flag (time infinity)
        bool it = false; // iteration limit flag
        if(endValue <= 0) // if the endValue specified by the user is non positive, then the infinite time flag is set
            ti = true;
        else
            cout << "endValue is " << endValue << endl; // otherwise, print the end value for confirmation to the user
        double epsilon = 0.05;
        bool checkEpsilon = false;
        bool checkRho = false;
        for(int i = 4; i < argc; i++) { // for all subsequent command line arguments (flags)
            string argument = argv[i];
            for (int j = 0; j < argument.length(); j++) { // set all strings fully to lower case format
                argument[j] = tolower(argument[j]);
            }
            // depending on the user input, set certain flags to true
            //if (argument =="-so" || argument == "-statesonly")
            //so = true;
            if (checkEpsilon) {
                epsilon = stod(argument);
                checkEpsilon = false;
            } else if (checkRho) {
                double rho = stod(argument);
                epsilon = (3.0 / (4.0*reactantsVector.size())) * (1.0 - sqrt((1.0+(rho/9.0)) / (1.0+rho)));
                checkEpsilon = false;
            } else {
                if (argument == "-fo" || argument == "-finalonly")
                    fo = true;
                else if (argument == "-it" || argument == "-useiter")
                    it = true;
                else if (argument == "-e" || argument == "-epsilon")
                    checkEpsilon = true;
                else if (argument == "-p" || argument == "-rho") {
                    checkRho = true;
                }
            }
        }

        if (epsilon > 1 || epsilon < 0) {
            cout << "An invalid epsilon or rho value was provided." << endl;
        } else {
            // Create directMethodSSA object and begin simulation
            boundedTauLeaping *btlAlgorithm = new boundedTauLeaping(moleculeAmounts, reactionRates, reactantsVector,
                                                                    stateChangeVector, endValue, fo, ti, it, epsilon);
            btlAlgorithm->start();

            // Get all necessary state from the directMethodSSA object
            vector<double> allTimes = btlAlgorithm->getAllTimes();
            vector<vector<int> > allStates = btlAlgorithm->getAllStates();
            vector<int> currentState = btlAlgorithm->getCurrentState();
            double currentTime = btlAlgorithm->getCurrentTime();
            int currentIteration = btlAlgorithm->getCurrentIteration();

            // Writing to output file depending on input flags
            ofstream outfile;
            outfile.open(outputFilePath);
            int iterationWidth = 15;
            char separator = ' ';
            if (fo) {
                //if(so) // final only and state only
                // outfile << "Final State:" << endl;
                //else { // final only with time
                if (it) {
                    outfile << "Final Iteration" << "\t\t\tFinal State" << endl;
                    outfile << left << setw(iterationWidth) << setfill(separator) << currentIteration << "\t\t\t";
                } else {
                    outfile << "Final Time" << "\t\t\tFinal State" << endl;
                    outfile << setprecision(5) << scientific << currentTime << "\t\t\t";
                }
                //}
                for (int i = 0; i < currentState.size(); i++) {
                    outfile << speciesList[i];
                    outfile << ": " << currentState[i];
                    outfile << "\t";
                }
            }
//        else if(so){ // states only (no time)
//            outfile << "Iteration" << "\t\t\tState" << endl;
//            for(int i = 0; i < allStates.size(); i++){
//                outfile << left << setw(iterationWidth) << setfill(separator) << i << "\t\t";
//                for(int j = 0; j < allStates[i].size(); j++){
//                    outfile << speciesList[j];
//                    outfile << ": " << allStates[i][j];
//                    outfile << "\t";
//                }
//                outfile << endl;
//            }
//        }
            else { // all states and all times
                outfile << "Iteration" << "\t\t\tTime(s)" << "\t\t\tState" << endl;
                for (int i = 0; i < allTimes.size(); i++) {
                    outfile << left << setw(iterationWidth) << setfill(separator) << i << "\t\t" << setprecision(5)
                            << scientific << allTimes[i] << "\t\t";
                    for (int j = 0; j < allStates[i].size(); j++) {
                        outfile << speciesList[j];
                        outfile << ": " << allStates[i][j];
                        outfile << "\t";
                    }
                    outfile << endl;
                }
            }
            outfile.close();
        }
    }
    else {
        cout << "There were issue(s) with your input file." << endl;
    }
    return 0;
}


