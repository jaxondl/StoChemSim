#ifndef BOUNDEDTAULEAPING_H
#define BOUNDEDTAULEAPING_H

#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include "../common/dependencyGraph.h"

class boundedTauLeaping {
private:
    dependencyGraph* dependency_graph;

    vector<double> reactionRates;
    vector<vector<pair<int, int> > > reactantsVector;
    vector<vector<pair<int, int> > > stateChangeVector;
    vector<int> currentState;
    vector<vector<int> > allStates;
    vector<double> allTimes;
    double currentTime;
    int currentIteration;
    double endValue;
    bool nonePossible;
    bool endInfinity;
    bool finalOnly;
    bool endByIteration;

    double epsilon;

    double getGammaRandomVariable(double a, double b);
    int getBinomialRandomVariable(int n, double p);

    vector<double> calculatePropensities();
    vector<int> calculateBounds(vector<double> propensities);
    vector<double> determineViolatingTimes(vector<int> bounds, vector<double> propensities);
    int determineFirstViolating(vector<double> violatingTimes);
    vector<int> determineReactionOccurrences(vector<int> bounds, vector<double> violatingTimes, int violatingIndex);

    void updateTime(double timeUntilNextReaction);
    void updateState(vector<int> reactionOccurrences);
    double getTotalPropensity();

public:
    boundedTauLeaping(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool finalOnly, bool endInfinity, bool endByIteration, double epsilon);

    vector<vector<int> > getAllStates();
    vector<double> getAllTimes();
    vector<int> getCurrentState();
    double getCurrentTime();
    int getCurrentIteration();

    void start();
};

#endif
