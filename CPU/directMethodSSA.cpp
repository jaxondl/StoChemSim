#include "directMethodSSA.h"

using namespace std;

directMethodSSA::directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double tEnd, bool statesOnly, bool finalOnly, bool tInfinity){
    this->reactionTree = new class reactionTree(moleculeAmounts, reactionRates, reactantsVector);
    this->dependencyGraph = new class dependencyGraph(stateChangeVector, reactantsVector);
    this->allStates.push_back(moleculeAmounts);
    this->allTimes.push_back(0);
    this->reactionRates = reactionRates; // k reaction constants
    this->reactantsVector = reactantsVector;
    this->stateChangeVector = stateChangeVector;
    this->currentState = moleculeAmounts;
    this->currentTime = 0;
    this->tEnd = tEnd;
    this->statesOnly = statesOnly;
    this->finalOnly = finalOnly;
    this->tInfinity = tInfinity;
    if (this->statesOnly) {
        this->tInfinity = true; 
    }
}

double directMethodSSA::getUniformRandomVariable(){
    random_device rd; // random seed
    mt19937 gen(rd()); // mersenne twister engine
    uniform_real_distribution<> dis(0.0, 1.0);
    double RV = dis(gen);
    return RV;
}

void directMethodSSA::updateTime(double timeUntilNextReaction){
    currentTime += timeUntilNextReaction;
}

void directMethodSSA::updateState(vector<vector<pair<int, int>>> stateChangeVector, int reactionIndex) {
    vector<pair<int, int>> chosenReactionChange = stateChangeVector[reactionIndex];
    for(pair<int, int> p: chosenReactionChange){
        currentState[p.first] += p.second; //update the state change at the associated index
    }

}

double directMethodSSA::getTimeUntilNextReaction(double propensity) {
    random_device rd;
    mt19937 gen(rd());
    exponential_distribution<> dis(propensity); //exponential approximation of propensity
    double RV = dis(gen);
    return RV;
}

double directMethodSSA::getTotalPropensity(){
    reactionTree::reactionNode root = reactionTree->reactionTreeArray[0];
    double totalPropensity = root.propensity + root.leftSum + root.rightSum;
    return totalPropensity;
}

vector<vector<int>> directMethodSSA::getAllStates(){
    return allStates;
}

vector<double> directMethodSSA::getAllTimes(){
    return allTimes;
}

vector<int> directMethodSSA::getCurrentState(){
    return currentState;
}

double directMethodSSA::getCurrentTime(){
    return currentTime;
}

void directMethodSSA::start(){
    while (getTotalPropensity() > 0.001 && (currentTime < tEnd || tInfinity)){
        if (!statesOnly) {
            double timeUntilNextReaction = getTimeUntilNextReaction(getTotalPropensity());
            if(currentTime + timeUntilNextReaction > tEnd && !tInfinity)
                break;
            updateTime(timeUntilNextReaction);
            allTimes.push_back(currentTime);
        }
        double uniformRV = getUniformRandomVariable();
        int reactionIndex = reactionTree->searchForNode(uniformRV);
        updateState(stateChangeVector, reactionIndex);
        if (!finalOnly) {
            allStates.push_back(currentState);
        }
        vector<int> dependentReactionIndices = dependencyGraph->getDependentReactions(reactionIndex);
        for(int reaction: dependentReactionIndices){
            reactionTree->updatePropensity(reaction, reactionRates[reaction], currentState, reactantsVector[reaction]);
        }
    }
    // If we pass tEnd, remove the last time and state before sending
}