#include "directMethodSSA.h"

using namespace std;

directMethodSSA::directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool statesOnly, bool finalOnly, bool endInfinity, bool endByIteration){
    this->reaction_tree = new class reactionTree(moleculeAmounts, reactionRates, reactantsVector); // create reaction tree
    this->dependency_graph = new class dependencyGraph(stateChangeVector, reactantsVector, moleculeAmounts); // create dependency graph
    this->allStates.push_back(moleculeAmounts); // stored initial state in the allStates vector
    this->allTimes.push_back(0); // store the initial time (0) in the allTimes vector
    this->reactionRates = reactionRates; // k reaction constants
    this->reactantsVector = reactantsVector;
    this->stateChangeVector = stateChangeVector;

    this->currentState = moleculeAmounts;
    this->currentTime = 0; // want to set the current time to 0 before beginning
    this->currentIteration = 0; // likewise set current iteration to 0 before beginning

    this->endValue = endValue; // all related booleans/flags are defined by the user
    this->statesOnly = statesOnly;
    this->finalOnly = finalOnly;
    this->endInfinity = endInfinity;
    this->endByIteration = endByIteration;

    random_device rd_uni; // random seed
    mt19937 gen_uni(rd_uni()); // mersenne twister engine
    this->gen_uni = gen_uni;
    random_device rd_exp; // random seed
    mt19937 gen_exp(rd_exp()); 
    this->gen_exp = gen_exp; // mersenne twister engine
    uniform_real_distribution<double> dis(0.0, 1);
    
    if (this->statesOnly && !this->endByIteration) { // if the user requests to only calculate the states but gives a positive finite value for time (i.e. did not set endByIteration to true) then the end time is automatically set to infinity (as finite time will be unapplicable)
        this->endInfinity = true;
    }
    //asdasdasdasdasdas
}

double directMethodSSA::getUniformRandomVariable() {
    //double totalPropensity = getTotalPropensity();
    uniform_real_distribution<double> dis(0.0, 1);
    double RV = dis(gen_uni); // obtain sampled value
    return RV;
}

void directMethodSSA::updateTime(double timeUntilNextReaction){
    currentTime += timeUntilNextReaction; // update time by adding the passed in time until next reaction to the existing time
}

void directMethodSSA::updateState(vector<vector<pair<int, int> > > stateChangeVector, int reactionIndex) {
    vector<pair<int, int> > chosenReactionChange = stateChangeVector[reactionIndex]; // given the chosen reaction, obtain all the changes to all the species
    for(pair<int, int> p: chosenReactionChange){
        currentState[p.first] += p.second; // update the state change by individually updating the amounts of the affected species
    }
}

double directMethodSSA::getTimeUntilNextReaction(double propensity) {
    exponential_distribution<> dis(propensity); // utilize exponential approximation of propensity
    double RV = dis(gen_exp); // obtain sampled value
    return RV;
}

double directMethodSSA::getTotalPropensity(){
    reactionTree::reactionNode root = reaction_tree->reactionTreeArray[0]; // propensity sum = root node propensity + rightSum + leftSum
    double totalPropensity = root.propensity + root.leftSum + root.rightSum;
    return totalPropensity;
}

vector<vector<int> > directMethodSSA::getAllStates(){return allStates;}

vector<double> directMethodSSA::getAllTimes(){return allTimes;}

vector<int> directMethodSSA::getCurrentState(){return currentState;}

double directMethodSSA::getCurrentTime(){return currentTime;}

int directMethodSSA::getCurrentIteration(){return currentIteration;}

void directMethodSSA::start(){
    // continue the simulation while the total propensity > 0 AND (the endInfinity flag is true OR the current time/iteration hasn't exceeded the inputted limit)
    while (getTotalPropensity() > 0 && ((!endByIteration && (currentTime < endValue || endInfinity)) || (endByIteration && (currentIteration < endValue || endInfinity)))){
        if (!statesOnly) { // only calculate if the user wants to also calculate the times (default)
            double timeUntilNextReaction = getTimeUntilNextReaction(getTotalPropensity()); // obtain the time until the next reaction
            if(currentTime + timeUntilNextReaction > endValue && !endInfinity) // if updating the time violates the finite end time value, terminate the simulation
                break;
            updateTime(timeUntilNextReaction); // otherwise, update the time
            allTimes.push_back(currentTime); // recorded the updated time
        }
        double RV = getUniformRandomVariable(); // obtain sampled value
        int reactionIndex = reaction_tree->searchForNode(RV); // search for reaction
        if(reactionIndex == -1){
            allTimes.pop_back();
            break;
        }
        updateState(stateChangeVector, reactionIndex); // update state/configuration
        if (!finalOnly) { // only save the state vectors of the iteration if the finalOnly flag is false
            allStates.push_back(currentState);
        }
        vector<int> dependentReactionIndices = dependency_graph->getDependentReactions(reactionIndex); // get affected reactions from dependency graph
        for(int reaction: dependentReactionIndices){
            reaction_tree->updatePropensity(reaction, reactionRates[reaction], currentState, reactantsVector[reaction]); // update propensity for every affected reaction
        }
        currentIteration++; // update iteration
    }
    cout << endl << "Final Iteration: " << currentIteration << endl;
    cout << "Final Time: " << currentTime << endl;
    cout << "Final Total Propensity: " << getTotalPropensity() << endl;
}
