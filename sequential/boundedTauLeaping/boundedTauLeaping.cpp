#include "boundedTauLeaping.h"

using namespace std;

boundedTauLeaping::boundedTauLeaping(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double endValue, bool finalOnly, bool endInfinity, bool endByIteration, double epsilon){
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
    this->finalOnly = finalOnly;
    this->endInfinity = endInfinity;
    this->endByIteration = endByIteration;
    this->epsilon = epsilon;
}

vector<double> boundedTauLeaping::calculatePropensities(){

}

vector<int> boundedTauLeaping::calculateBounds(vector<double> propensities){

}

vector<double> boundedTauLeaping::determineViolatingTimes(vector<int> bounds, vector<double> propensities){

}

int boundedTauLeaping::determineFirstViolating(vector<double> violatingTimes){

}

vector<int> boundedTauLeaping::determineReactionOccurrences(vector<int> bounds, vector<double> violatingTimes, int violatingIndex){

}

void boundedTauLeaping::updateTime(double timeUntilNextReaction){

}

void boundedTauLeaping::updateState(vector<vector<pair<int, int>>> stateChangeVector, int reactionIndex){

}
double boundedTauLeaping::getTotalPropensity(){
    
}

vector<vector<int>> boundedTauLeaping::getAllStates(){return allStates;}

vector<double> boundedTauLeaping::getAllTimes(){return allTimes;}

vector<int> boundedTauLeaping::getCurrentState(){return currentState;}

double boundedTauLeaping::getCurrentTime(){return currentTime;}

int boundedTauLeaping::getCurrentIteration(){return currentIteration;}
    
void boundedTauLeaping::start(){
    
}