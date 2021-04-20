#ifndef BOUNDEDTAULEAPING_H
#define BOUNDEDTAULEAPING_H

using namespace std;

#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <algorithm>
#include <climits>

class boundedTauLeaping {
private:
    vector<vector<int> > allStates;
    vector<double> allTimes;
    vector<double> reactionRates;
    vector<vector<pair<int, int> > > reactantsVector;
    vector<vector<pair<int, int> > > stateChangeVector;

    vector<int> currentState;
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
    vector<int> calculateBounds();
    vector<double> determineViolatingTimes(vector<int> firingBounds, vector<double> propensities);
    int determineFirstViolating(vector<double> violatingTimes);
    vector<int> determineReactionOccurrences(vector<int> firingBounds, vector<double> violatingTimes, int violatingIndex);

    void updateTime(double tau);
    void updateState(vector<int> reactionOccurrences);

public:
    boundedTauLeaping(vector<int> initialState, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool finalOnly, bool endInfinity, bool endByIteration, double epsilon);

    vector<vector<int> > getAllStates();
    vector<double> getAllTimes();
    vector<int> getCurrentState();
    double getCurrentTime();
    int getCurrentIteration();

    void start();
};

#endif
