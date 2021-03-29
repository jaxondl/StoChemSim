#include "boundedTauLeaping.h"
#include <cmath> // for absolute value
#include <random>

using namespace std;

boundedTauLeaping::boundedTauLeaping(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool finalOnly, bool endInfinity, bool endByIteration, double epsilon){
    this->dependency_graph = new class dependencyGraph(stateChangeVector, reactantsVector, moleculeAmounts);
    this->allStates.push_back(moleculeAmounts);
    this->allTimes.push_back(0);
    this->reactionRates = reactionRates; // k reaction constants
    this->reactantsVector = reactantsVector;
    this->stateChangeVector = stateChangeVector;
    this->currentState = moleculeAmounts;
    this->currentTime = 0;
    this->currentIteration = 0;
    this->endValue = endValue;
    this->nonePossible = false;
    this->finalOnly = finalOnly;
    this->endInfinity = endInfinity;
    this->endByIteration = endByIteration;
    this->epsilon = epsilon;
}

//a is alpha, b is beta; both must be positive
double boundedTauLeaping::getGammaRandomVariable(double a, double b){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    gamma_distribution<double> distribution(a, b);
    double g = distribution(generator);
    // cout << "gamma RV: " << g << endl;
    return g;
}

int boundedTauLeaping::getBinomialRandomVariable(int n, double p){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    binomial_distribution<double> distribution(n, p);
    double b = distribution(generator);
    // cout << "binomial RV: " << b << endl;
    return b;
}

vector<double> boundedTauLeaping::calculatePropensities(){
    // Propesnity calculation for each reaction = k * (amount(reactantMolecule1)) * (amount(reactantMolecule1)-1)...(amount(reactantMolecule1)-reactantCoefficient + 1) * ...
    nonePossible = true;
    vector<double> propensities;
    for (int i = 0; i < reactionRates.size(); i++) {
        double propensity = reactionRates[i];
        for (pair<int, int> reactant : reactantsVector[i]) {
            for (int j = 0; j < reactant.second; j++) {
                propensity *= (currentState[reactant.first] - i);
            }
        }
        if (propensity > 0) {
            nonePossible = false;
        }
        propensities.push_back(propensity);
    }
    return propensities;
}

vector<int> boundedTauLeaping::calculateBounds(vector<double> propensities){
    vector<int> result;
    for (int i = 0; i < propensities.size(); i++) {
        vector<pair<int, int> > currentStateChangeVector = stateChangeVector[i]; // looks through each reaction's state changes

        // |bj * vij| > eij * xi, hence eij = 0.25, xi = 2, then eij * xi = 0
        // bj > |eij * xi / vij|
        // bj = ceil(|eij * xi / vij|)

        int bound = INT_MAX;
        for (pair<int, int> molecule : currentStateChangeVector) {
            int eTimesX = epsilon * currentState[molecule.first];
            int current = ceil(abs(eTimesX / molecule.second));
            
            if (current < bound) {
                bound = current;
            }
        }
        result.push_back(bound);
    }
    return result;
}

vector<double> boundedTauLeaping::determineViolatingTimes(vector<int> bounds, vector<double> propensities){
    vector<double> result;
    for (int i = 0; i < bounds.size(); i++) {
        result.push_back(getGammaRandomVariable(bounds[i], propensities[i]));
    }
    return result;
}

int boundedTauLeaping::determineFirstViolating(vector<double> violatingTimes){
    vector<double>::iterator iter = min_element(violatingTimes.begin(), violatingTimes.end());
    return distance(violatingTimes.begin(), iter);
}

vector<int> boundedTauLeaping::determineReactionOccurrences(vector<int> bounds, vector<double> violatingTimes, int violatingIndex) {
    vector<int> result;
    double firstViolatingTime = violatingTimes[violatingIndex];
    for (int i = 0; i < violatingTimes.size(); i++) {
        int nj;
        if (i == violatingIndex) {
            nj = bounds[violatingIndex];
        }
        else {
            double p = firstViolatingTime/violatingTimes[i];
            nj = getBinomialRandomVariable(bounds[i] - 1, p);
        }
        result.push_back(nj);
    }
    return result;
}

void boundedTauLeaping::updateTime(double timeUntilNextReaction){
    currentTime += timeUntilNextReaction;
}

//Effect the leap by replacing ~x[j] with ~x[j] + ~ν[j]*n[j]
//notation: ~x means the x vector (v is the state change vector, x is the current state vector)
//nj is a binomial RV from with parameters (bj - 1, τ / τj) (aka the value returned by determineReactionOccurrences())
//b is the bounds vector (returned by calculateBounds())
void boundedTauLeaping::updateState(vector<int> reactionOccurrences) {
    for (int i = 0; i < stateChangeVector.size(); i++) {
        vector<pair<int, int> > chosenReactionChange = stateChangeVector[i]; // Vj vector
        int numberOccurrences = reactionOccurrences[i]; // n
        for (pair<int, int> p: chosenReactionChange) {
            currentState[p.first] += p.second*numberOccurrences; // x[j] += Vj*Nj
            // update the state change by individually updating the amounts of the affected species, for each time the reaction occurred
        }
    }
}
double boundedTauLeaping::getTotalPropensity(){

}

vector<vector<int> > boundedTauLeaping::getAllStates(){return allStates;}

vector<double> boundedTauLeaping::getAllTimes(){return allTimes;}

vector<int> boundedTauLeaping::getCurrentState(){return currentState;}

double boundedTauLeaping::getCurrentTime(){return currentTime;}

int boundedTauLeaping::getCurrentIteration(){return currentIteration;}
    
void boundedTauLeaping::start(){
    vector<double> props = calculatePropensities(); // should determine nonePossible
    vector<int> firingBounds;
    vector<double> violatingTimes;
    int firstViolatingIndex;
    vector<int> reactionOccurrences;

    while (!nonePossible && ((!endByIteration && (currentTime < endValue || endInfinity)) || (endByIteration && (currentIteration < endValue || endInfinity)))){
        firingBounds = calculateBounds(props); // Step 1b
        violatingTimes = determineViolatingTimes(firingBounds, props); // Step 2
        firstViolatingIndex = determineFirstViolating(violatingTimes); // Step 3
        reactionOccurrences = determineReactionOccurrences(firingBounds, violatingTimes, firstViolatingIndex); // Step 4
        
        if(currentTime + violatingTimes[firstViolatingIndex] > endValue && !endInfinity) // if updating the time violates the finite end time value, terminate the simulation
            break;
        updateTime(violatingTimes[firstViolatingIndex]); // Step 5 effecting the leap: time update
        updateState(reactionOccurrences); // Step 5 effecting the leap : state molevule amounts
         // recorded the updated time
        if (!finalOnly) { // only save the state and time vectors of the iteration if the finalOnly flag is false
            allTimes.push_back(currentTime);
            allStates.push_back(currentState);
        }
        
        currentIteration++; // update iteration
        props = calculatePropensities(); // actually step 1a, but we caclulated outside the while loop
    }    
}
