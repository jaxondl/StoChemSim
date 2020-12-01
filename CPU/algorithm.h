#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <iostream>
#include <random>
#include <vector>
#include "dependency.h"
#include "tree.h"


class algorithm {
private:
    tree* reactionTree; 
    dependency* dependencyGraph;

    vector<double> reactionRates;
    vector<vector<pair<int, int>>> reactantsVector;
    vector<vector<pair<int, int>>> stateChangeVector;
    vector<int> currentState;
    double currentTime;
    double t_end;

    //functions
    double getUniformRandomVariable();
    double getTimeUntilNextReaction(double propensity); //returns Exp(propensity)
    void updateTime(double timeUntilNextReaction);
    void updateState(vector<vector<pair<int, int>>> stateChangeVector, int reactionIndex);
    double getTotalPropensity();

public:
    algorithm(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double t_end);
    void start();
};

#endif //ALGORITHM_H
