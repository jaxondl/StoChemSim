#ifndef DIRECTMETHODSSA_H
#define DIRECTMETHODSSA_H

#include <iostream>
#include <random>
#include <vector>
#include "dependencyGraph.h"
#include "reactionTree.h"


class directMethodSSA {
private:
    reactionTree* reaction_tree; // reaction tree necessary for this simulation algorithm
    dependencyGraph* dependency_graph; // dependency graph necessary for this simulation algorithm

    vector<double> reactionRates; // reaction rates (k) for each reaction
    vector<vector<pair<int, int> > > reactantsVector; // reactants for each reaction
    vector<vector<pair<int, int> > > stateChangeVector; // state change vectors for each reaction
    vector<int> currentState; // the current state / amounts of species in CRN
    vector<vector<int> > allStates; // vector containing all states calculated by the simulation algorithm
    vector<double> allTimes; // vector containing the time it has taken to complete each iteration and all preceding iterations

    double currentTime; // current time tracked for the simulation
    int currentIteration; // current iteration tracked fo the simulation
    double endValue; // end value will either be a ending time or and ending iteration, depending on whether the user sets the endByIteration flag to true. The simulation will stop before it crosses the endValue

    bool endInfinity; // if set to true, the simulation will continue until no more reactions can take place
    bool statesOnly; // if set to true, the simulation will forego any calculations of times (as well as any updating the time) and only prin out the states 
    bool finalOnly; // if set to true, the simulation will only output/print the final iteration of the simulation
    bool endByIteration; // if set to true, the endValue will be used to determine the ending iteration instead of the ending time

    double getUniformRandomVariable(); // obtain a uniform RV dist sample value
    double getTimeUntilNextReaction(double propensity); // obtain the time until the next reaction (in order to update the time)
    void updateTime(double timeUntilNextReaction); // update the existing time
    void updateState(vector<vector<pair<int, int> > > stateChangeVector, int reactionIndex); // update the state of the species
    double getTotalPropensity(); // obtain the total propensity across all reactions

public:
    directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool statesOnly, bool finalOnly, bool endInfinity, bool endByIteration); // constructor
    vector<vector<int> > getAllStates(); // obtain vector of all states
    vector<double> getAllTimes(); // obtain vector of all times
    vector<int> getCurrentState(); // obtain vector of the end/most recent state
    double getCurrentTime(); // obtain end/most recent time
    int getCurrentIteration(); // obtain end/most recent iteration
    void start(); // begin CRN SSA simulation
};

#endif
