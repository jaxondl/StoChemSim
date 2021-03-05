#ifndef DIRECTMETHODSSA_H
#define DIRECTMETHODSSA_H

#include <iostream>
#include <random>
#include <vector>
#include "dependencyGraph.h"
#include "reactionTree.h"


class directMethodSSA {
private:
    reactionTree* reaction_tree; 
    dependencyGraph* dependency_graph;

    vector<double> reactionRates;
    vector<vector<pair<int, int>>> reactantsVector;
    vector<vector<pair<int, int>>> stateChangeVector;
    vector<int> currentState;
    vector<vector<int>> allStates;
    vector<double> allTimes;
    double currentTime;
    int currentIteration;
    double endValue;
    bool endInfinity;
    bool statesOnly;
    bool finalOnly;
    bool endByIteration;

    double getUniformRandomVariable();
    double getTimeUntilNextReaction(double propensity);
    void updateTime(double timeUntilNextReaction);
    void updateState(vector<vector<pair<int, int>>> stateChangeVector, int reactionIndex);
    double getTotalPropensity();

public:
    directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double endValue, bool statesOnly, bool finalOnly, bool endInfinity, bool endByIteration);
    vector<vector<int>> getAllStates();
    vector<double> getAllTimes();
    vector<int> getCurrentState();
    double getCurrentTime();
    int getCurrentIteration();
    void start();
};

#endif
