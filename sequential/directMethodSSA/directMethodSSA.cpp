#include "directMethodSSA.h"

using namespace std;

directMethodSSA::directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double endValue, bool statesOnly, bool finalOnly, bool endInfinity, bool endByIteration){
    this->reaction_tree = new class reactionTree(moleculeAmounts, reactionRates, reactantsVector);
    this->dependency_graph = new class dependencyGraph(stateChangeVector, reactantsVector);
    this->allStates.push_back(moleculeAmounts);
    this->allTimes.push_back(0);
    this->reactionRates = reactionRates; // k reaction constants
    this->reactantsVector = reactantsVector;
    this->stateChangeVector = stateChangeVector;
    this->currentState = moleculeAmounts;
    this->currentTime = 0;
    this->currentIteration = 0;
    this->endValue = endValue;
    this->statesOnly = statesOnly;
    this->finalOnly = finalOnly;
    this->endInfinity = endInfinity;
    this->endByIteration = endByIteration;
    if (this->statesOnly && !this->endByIteration) {
        this->endInfinity = true;
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
    reactionTree::reactionNode root = reaction_tree->reactionTreeArray[0];
    double totalPropensity = root.propensity + root.leftSum + root.rightSum;
    return totalPropensity;
}

vector<vector<int>> directMethodSSA::getAllStates(){return allStates;}

vector<double> directMethodSSA::getAllTimes(){return allTimes;}

vector<int> directMethodSSA::getCurrentState(){return currentState;}

double directMethodSSA::getCurrentTime(){return currentTime;}

int directMethodSSA::getCurrentIteration(){return currentIteration;}

void directMethodSSA::start(){
    while (getTotalPropensity() > 0.001 && ((!endByIteration && (currentTime < endValue || endInfinity)) || (endByIteration && (currentIteration < endValue || endInfinity)))){
        if (!statesOnly) {
            double timeUntilNextReaction = getTimeUntilNextReaction(getTotalPropensity());
            if(currentTime + timeUntilNextReaction > endValue && !endInfinity)
                break;
            updateTime(timeUntilNextReaction);
            allTimes.push_back(currentTime);
        }
        double uniformRV = getUniformRandomVariable();
        int reactionIndex = reaction_tree->searchForNode(uniformRV);
        updateState(stateChangeVector, reactionIndex);
        if (!finalOnly) {
            allStates.push_back(currentState);
        }
        vector<int> dependentReactionIndices = dependency_graph->getDependentReactions(reactionIndex);
        for(int reaction: dependentReactionIndices){
            reaction_tree->updatePropensity(reaction, reactionRates[reaction], currentState, reactantsVector[reaction]);
        }
        currentIteration++;
    }
    // If we pass endValue, remove the last time and state before sending
}