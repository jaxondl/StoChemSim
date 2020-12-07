#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include "inputVerifier.h"


//#include <bits/stdc++.h> 
using namespace std;

struct Reaction {
    vector<string> reactantTypes;
    vector<string> productTypes;
    vector<int> reactantCoefficients;
    vector<int> productCoefficients;
    int k; //reaction constant
    bool usable = true;
    long propensity;
};

struct State {
    unordered_map<string, int> moleculeCoefficientMap; 
};

void SSA(State currentState, vector<Reaction> reactions) {
    cout << "SSA begin" << endl;
    bool endSimulationFlag = false;
    int iterations = 0;
    while (true) {
        // Check if the simulation is over (no reaction can take place)
        vector<Reaction> validReactions;
        int propensitySum = 0;
        for (auto reaction : reactions) {
            bool reactionUsable = true;
            int coefficientIndex = 0;
            int propensityProduct = reaction.k;
            for(auto itr : reaction.reactantTypes) {
                int moleculeAmount = currentState.moleculeCoefficientMap[itr];
                propensityProduct *= moleculeAmount;
                int reactantCoef = reaction.reactantCoefficients.at(coefficientIndex);
                if (reactantCoef > moleculeAmount) {
                    reactionUsable = false;
                }
                coefficientIndex++;
            }
            reaction.usable = reactionUsable; // sets whether the reaction is usable or not at every iteration
            if(reactionUsable) {
                reaction.propensity = propensityProduct;
                propensitySum += propensityProduct;
                validReactions.push_back(reaction);
            }
        }
        if (validReactions.empty()){
            break;
        }
        vector<double> probabilityEndpoints;
        double probabilityPoint = 0;
        for (auto reaction:validReactions) { // map out number line
            double propensity = reaction.propensity;
            double probability = propensity / propensitySum;
            probabilityPoint += probability;
            probabilityEndpoints.push_back(probabilityPoint);
        }

        random_device rd; // random seed
        mt19937 gen(rd()); // mersenne twister engine
        uniform_real_distribution<> dis(0.0, 1.0);
        double RV = dis(gen);
        cout << "Chosen U: " << RV << endl;
        int chosenReactionIndex;
        for(int i = 0; i < probabilityEndpoints.size(); i++){
            if(RV <= probabilityEndpoints[i]){
                chosenReactionIndex = i;
                break;
            }
        }
        Reaction chosenReaction = validReactions.at(chosenReactionIndex);
        int coefficientIndex = 0;
        for(auto itr : chosenReaction.reactantTypes) {
            int reactantCoef = chosenReaction.reactantCoefficients.at(coefficientIndex);
            currentState.moleculeCoefficientMap[itr] -= reactantCoef; // removes the relevant reactant from the chosen reaction from the state
            coefficientIndex++;
        }
        coefficientIndex = 0;
        for(auto itr : chosenReaction.productTypes) {
            int productCoef = chosenReaction.productCoefficients.at(coefficientIndex);
            currentState.moleculeCoefficientMap[itr] += productCoef; // adds the relevant product from the chosen reaction from the state
            coefficientIndex++;
        }

        for (auto molecule : currentState.moleculeCoefficientMap) {
            cout << molecule.first << " " << molecule.second << endl;
        }
        validReactions.clear();
        probabilityEndpoints.clear();
        iterations++;
    }
}

int main() {
    inputVerifier *iv = new inputVerifier;
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\sample_input_SSA_file.txt");
    bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\faulty_sample_input_SSA_file_1.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\faulty_sample_input_SSA_file_2.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\faulty_sample_input_SSA_file_3.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\faulty_sample_input_SSA_file_4.txt");
    //bool safeToRun = iv->verifyFile("C:\\Users\\ccort\\CLionProjects\\ssa-implementation\\faulty_sample_input_SSA_file_5.txt");

    if (safeToRun) {
        State initialState;
        initialState.moleculeCoefficientMap["A"] = 5;
        initialState.moleculeCoefficientMap["B"] = 3;
        initialState.moleculeCoefficientMap["C"] = 0;

        Reaction reaction1;
        reaction1.reactantTypes.push_back("A");
        reaction1.productTypes.push_back("B");
        reaction1.reactantCoefficients.push_back(1);
        reaction1.productCoefficients.push_back(1);
        reaction1.k = 2;

        Reaction reaction2;
        reaction2.reactantTypes.push_back("A");
        reaction2.reactantTypes.push_back("B");
        reaction2.productTypes.push_back("C");
        reaction2.reactantCoefficients.push_back(1);
        reaction2.reactantCoefficients.push_back(1);
        reaction2.productCoefficients.push_back(1);
        reaction2.k = 3;

        vector<Reaction> reactions;
        reactions.push_back(reaction1);
        reactions.push_back(reaction2);

        SSA(initialState, reactions);
    }
    return 0;
}