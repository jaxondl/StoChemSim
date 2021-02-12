#ifndef DIRECTMETHODSSA_H
#define DIRECTMETHODSSA_H

#include <iostream>
#include <random>
#include <vector>
#include "dependencyGraph.h"
#include "reactionTree.h"


class directMethodSSA {
private:
    reactionTree* reactionTree; 
    dependencyGraph* dependencyGraph;

    vector<double> reactionRates;
    vector<vector<pair<int, int>>> reactantsVector;
    vector<vector<pair<int, int>>> stateChangeVector;
    vector<int> currentState;
    vector<vector<int>> allStates;
    vector<double> allTimes;
    double currentTime;
    double tEnd;
    bool tInfinity;
    bool statesOnly;
    bool finalOnly;

    double getUniformRandomVariable();
    double getTimeUntilNextReaction(double propensity);
    void updateTime(double timeUntilNextReaction);
    void updateState(vector<vector<pair<int, int>>> stateChangeVector, int reactionIndex);
    double getTotalPropensity();

public:
    directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double t_end, bool statesOnly, bool finalOnly, bool tInfinity);
    vector<vector<int>> getAllStates();
    vector<double> getAllTimes();
    vector<int> getCurrentState();
    double getCurrentTime();
    void start();
};

#endif
