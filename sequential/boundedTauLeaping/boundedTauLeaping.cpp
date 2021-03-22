#include "boundedTauLeaping.h"
#include <cmath> // for absolute value

using namespace std;

boundedTauLeaping::boundedTauLeaping(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool finalOnly, bool endInfinity, bool endByIteration, double epsilon){
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

//a is alpha, b is beta; both must be positive
double boundedTauLeaping::getGammaRandomVariable(double a, double b){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::gamma_distribution<double> distribution(a, b);
    double g = distribution(generator);
    cout << "gamma RV: " << g << endl;
    return g;
}

double boundedTauLeaping::getBinomialRandomVariable(int n, double p){
}

vector<double> boundedTauLeaping::calculatePropensities(){

}

vector<int> boundedTauLeaping::calculateBounds(vector<double> propensities){
    vector<int> result;
    for(int i = 0; i < propensities.size(); i++){
        int bound = 0; // set the integer to the minimum non-negative value
        bool withinBoundary = true;
        vector<pair<int, int>> currentStateChangeVector = stateChangeVector[i]; // looks through each reaction's state changes
        while(withinBoundary){
            bound += 1;
            for(auto entry : currentStateChangeVector){ // entry is Vij in the paper : the first of the pair is the molecule id, the second is its quantity
                if(abs(entry.second * bound) > (epsilon * currentState[entry.first])){ // if we reach a violation, we stop increasing the bound and stop the loop
                    withinBoundary = false;
                    break;
                }
            }
        }
        result.push_back(bound);
    }
    return result;
}

vector<double> boundedTauLeaping::determineViolatingTimes(vector<int> bounds, vector<double> propensities){
    vector<double> result;
    for(int i = 0; i < bounds.size(); i++){
        result.push_back(getGammaRandomVariable(bounds[i], propensities[i]));
    }
    return result;
}

int boundedTauLeaping::determineFirstViolating(vector<double> violatingTimes){
    vector<double>::iterator iter = min_element(violatingTimes.begin(), violatingTimes.end());
    return distance(violatingTimes.begin(), iter) - 1;
}

vector<int> boundedTauLeaping::determineReactionOccurrences(vector<int> bounds, vector<double> violatingTimes, int violatingIndex){

}

void boundedTauLeaping::updateTime(double timeUntilNextReaction){
    currentTime += timeUntilNextReaction;
}

void boundedTauLeaping::updateState(vector<vector<pair<int, int> > > stateChangeVector, int reactionIndex){

}
double boundedTauLeaping::getTotalPropensity(){

}

vector<vector<int> > boundedTauLeaping::getAllStates(){return allStates;}

vector<double> boundedTauLeaping::getAllTimes(){return allTimes;}

vector<int> boundedTauLeaping::getCurrentState(){return currentState;}

double boundedTauLeaping::getCurrentTime(){return currentTime;}

int boundedTauLeaping::getCurrentIteration(){return currentIteration;}
    
void boundedTauLeaping::start(){
    
}
