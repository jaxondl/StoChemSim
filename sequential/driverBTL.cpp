#include "common/decoder.h"
#include "common/inputVerifier.h"
#include "boundedTauLeaping/boundedTauLeaping.h"

#include <boost/program_options.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

using namespace std;
namespace po = boost::program_options;


int main(int argc, char** argv) {

   po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("endtime,t", po::value<double>()->default_value(0), "specify endtime if ending by time, nonpositive value ends at inifinty, do not set both this and enditer")
        ("enditer,i", po::value<int>(), "specify end iteration if ending by iteration, nonpositive value ends at inifinty, do not set both this and endtime")
        ("finalonly,f", po::value<bool>()->default_value(false), "set to true to only save final state and time")
        ("epsilon,e", po::value<double>()->default_value(0), "set epsilon value; default is calculated by default rho")
        ("rho,p", po::value<double>()->default_value(0.25), "set rho value; used to calculate epsilon; default is 0.25")
        ("out,o", po::value<string>(), "output file if specified")
        ("input-file", po::value<string>(), "input file, MUST be specified")
        ;
    po::positional_options_description p;
    p.add("input-file", 1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
//    po::store(po::parse_command_line(argc, argv, desc), vm);
//    po::notify(vm); 

    bool endByTime = false;
    bool endByIter = false;

    if (vm.count("input-file")) {
        cout << "Input file path is: " << vm["input-file"].as< string >() << "\n";
    }
    else {
        cout << "Input file path was not set.\n";
        return 1;
    }
    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    if (vm.count("endtime")) {
        cout << "Endtime: " << vm["endtime"].as< double >() << "\n";
        endByTime = true;
    }
    if (vm.count("enditer")) {
        cout << "End Iter: " << vm["enditer"].as< int >() << "\n";
        endByIter = true;
    }
    if (vm.count("finalonly")) {
        cout << "Final Only: " << vm["finalonly"].as< bool >() << "\n";
        endByIter = true;
    }
    if (vm.count("epsilon")) {
        cout << "Epsilon: " << vm["epsilon"].as< double >() << "\n";
    }
    if (vm.count("rho")) {
        cout << "Rho: " << vm["rho"].as< double >() << "\n";
    }
    
    inputVerifier *iv = new inputVerifier(); // create input verifier for validating input file
    string inputFilePath = argv[1]; // example: C:\\Users\\Isaac\\CLionProjects\\SeniorDesign\\crn-ssa-wolfram-pkg\\sequential\\inputs\\sample_input_SSA_file.txt
    //string inputFilePath = vm["input-file"].as< string >();
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
        double rho = 0.25;
        double epsilon = (3.0 / (4.0*reactantsVector.size())) * (1.0 - sqrt((1.0+(rho/9.0)) / (1.0+rho)));
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
                rho = stod(argument);
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
            cout << "Epsilon value is " << epsilon << endl;
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
