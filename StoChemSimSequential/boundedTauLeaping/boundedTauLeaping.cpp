#include "boundedTauLeaping.h"

using namespace std;

boundedTauLeaping::boundedTauLeaping(vector<int> initialState, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool finalOnly, bool endInfinity, bool endByIteration, double epsilon) {
    this->allStates.push_back(initialState);
    this->allTimes.push_back(0);
    this->reactionRates = reactionRates;
    this->reactantsVector = reactantsVector;
    this->stateChangeVector = stateChangeVector;
    this->currentState = initialState;
    this->currentTime = 0;
    this->currentIteration = 0;
    this->endValue = endValue;
    this->finalOnly = finalOnly;
    this->endInfinity = endInfinity;
    this->endByIteration = endByIteration;
    this->epsilon = epsilon;
}

// Draws from gamma distribution defined by a as alpha and b as beta, both must be positive
double boundedTauLeaping::getGammaRandomVariable(double a, double b) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    //if (1/b <= -0){    cout << "ERROR" << endl;}
    gamma_distribution<double> distribution(a, 1/b); // 1/b BECAUSE USE SCALE, NOT RATE
    double g = distribution(generator); // gets rid of potential -0 values and hence -inf gamma dist values
    //cout << "a: " << a << " b: " << b << " gamma RV: " << g << endl;
    return g;
}

// Draws from binomial distribution defined by n experiments (positive) and success probability p (between 0 and 1)
int boundedTauLeaping::getBinomialRandomVariable(int n, double p) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    binomial_distribution<int> distribution(n, p);
    int b = distribution(generator);
    //cout << "binomial RV: " << b << endl;
    return b;
}

// Propensity calculation for each reaction = k * (amount(reactantSpecies1)) * (amount(reactantSpecies1)-1)...(amount(reactantSpecies1)-reactantCoefficient + 1) * ...
vector<double> boundedTauLeaping::calculatePropensities() {
    nonePossible = true;
    vector<double> propensities;
    double total = 0;
    // compute propensity for every reaction
    for (int i = 0; i < reactionRates.size(); i++) {
        double propensity = reactionRates[i];
        // start with reaction rate and multiply by reactant amounts according to propensity formula
        for (pair<int, int> reactant : reactantsVector[i]) {
            for (int j = 0; j < reactant.second; j++) {
                if(currentState[reactant.first] < 0){
                    nonePossible = true;
                    cout << "Ran into negative state species counts!" << " Propensity: " << propensity << endl;
                    // throw new exception. Try catch block in the driver itself
                    return propensities;
                }

                if((currentState[reactant.first] - j) < 0) break;
                propensity *= (currentState[reactant.first] - j);
                //if(propensity == 0) {cout << "YO THE PROP IS " << propensity << endl;}
            }
        }
        // determine if all propensities are 0 (no reactions possible, end simulation)
        if (propensity > 0) {
            nonePossible = false;
        }

        // if (propensity < 0){
        //     nonePossible = true;
        //     cout << "Ran into negative propensity!" << " Propensity: " << propensity << endl;
        //     // throw new exception. Try catch block in the driver itself
        //     break;
        // }

        propensities.push_back(propensity);
        total += propensity;
        //cout << "PROP IS " << propensity << " ";
    }
    cout << "Total Propensity: " << total << endl << endl;
    return propensities;
}

// Firing bound for each reaction is minimum bj such that |bj * vij| > epsilon * xi for some species i
// This implies bj = ceil(|epsilon * xi / vij|) for each affected species
// For each reaction, minimum bj across all affected species is the firing bound
vector<int> boundedTauLeaping::calculateBounds() {
    vector<int> firingBounds;
    // compute firing bound for every reaction
    for (int i = 0; i < stateChangeVector.size(); i++) {
         // use state change vector to determine affected species
        vector<pair<int, int> > currentStateChangeVector = stateChangeVector[i];

        int minBound = INT_MAX;
        // for all species affected by reaction, determine bj and keep track of minimum
        for (pair<int, int> species : currentStateChangeVector) {
            int currentBound = ceil( abs( epsilon * currentState[species.first] / species.second));
            if(currentBound == 0){
                currentBound = 1;
            }
            //cout << "curBound is " << currentBound << "; ep*curState = " << epsilon * currentState[species.first] << "; speciesSecond change = " << species.second << endl;
            if (currentBound < minBound) {
                minBound = currentBound;
            }
        }
        firingBounds.push_back(minBound);
        //cout << "MIN BOUND IS " << minBound << " ";
    }
    //cout << endl;
    return firingBounds;
}

// Determines violating times Tj by drawing from gamma distribution for every reaction using firing bounds and propensities
vector<double> boundedTauLeaping::determineViolatingTimes(vector<int> firingBounds, vector<double> propensities) {
    vector<double> violatingTimes;
    // draw from gamma distribution for every reaction
    for (int i = 0; i < firingBounds.size(); i++) {
        violatingTimes.push_back(getGammaRandomVariable(firingBounds[i], propensities[i]));
    }
    return violatingTimes;
}

// Determines first violating time index via argmin
int boundedTauLeaping::determineFirstViolating(vector<double> violatingTimes) {
    vector<double>::iterator iter = min_element(violatingTimes.begin(), violatingTimes.end());
    return distance(violatingTimes.begin(), iter);
}

// Determines number of occurrences nj for each reaction during leap T
// For violating reaction, nj = bj
// For all other reactions, nj is drawn from binomial distribution with n = bj - 1 and p = T / Tj
vector<int> boundedTauLeaping::determineReactionOccurrences(vector<int> firingBounds, vector<double> violatingTimes, int violatingIndex) {
    vector<int> result;
    double firstViolatingTime = violatingTimes[violatingIndex];
    // determine number of occurrences for every reaction during leap
    for (int i = 0; i < violatingTimes.size(); i++) {
        int nj;
        // if violating reaction, nj = bj
        if (i == violatingIndex) {
            nj = firingBounds[violatingIndex];
        }
        // else, draw from binomial distribution to determine occurrences
        else {
            double p = firstViolatingTime/violatingTimes[i];
            nj = getBinomialRandomVariable(firingBounds[i] - 1, p);
            //cout << "firing bounds: " << firingBounds[i] << " p is " << p << " nj is " << nj << endl;
        }
        result.push_back(nj);
        //cout << "change " << nj << " times ";
    }
    //cout << endl;
    return result;
}

// Updates time by leap amount tau
void boundedTauLeaping::updateTime(double tau) {
    currentTime += tau;
}

// Effects the leap by replacing ~x[j] with ~x[j] + ~ν[j] * n[j]
// Notation: ~x means the x vector (v is the state change vector, x is the current state vector)
void boundedTauLeaping::updateState(vector<int> reactionOccurrences) {
    // update state for each reaction based on number of occurrences
    for (int i = 0; i < reactionOccurrences.size(); i++) {
        vector<pair<int, int> > chosenReactionChange = stateChangeVector[i]; // ~v[j] vector
        int numberOccurrences = reactionOccurrences[i]; // n[j]
        // update the state change by individually updating the amounts of the affected species for each time the reaction occurred
        for (pair<int, int> species: chosenReactionChange) {
            currentState[species.first] += species.second*numberOccurrences;
        }
    }
}

vector<vector<int> > boundedTauLeaping::getAllStates() {return allStates;}

vector<double> boundedTauLeaping::getAllTimes() {return allTimes;}

vector<int> boundedTauLeaping::getCurrentState() {return currentState;}

double boundedTauLeaping::getCurrentTime() {return currentTime;}

int boundedTauLeaping::getCurrentIteration() {return currentIteration;}
    
void boundedTauLeaping::start() {
    vector<double> props = calculatePropensities(); // Step 1a, sets nonePossible
    vector<int> firingBounds;
    vector<double> violatingTimes;
    int firstViolatingIndex;
    vector<int> reactionOccurrences;

    // run until no reactions are possible or time/iteration termination condition is met
    while (!nonePossible && ((!endByIteration && (currentTime < endValue || endInfinity)) || (endByIteration && (currentIteration < endValue || endInfinity)))) {
        cout << "Leap Iteration: " << currentIteration << endl;
        firingBounds = calculateBounds(); // Step 1b
        violatingTimes = determineViolatingTimes(firingBounds, props); // Step 2
        firstViolatingIndex = determineFirstViolating(violatingTimes); // Step 3
        reactionOccurrences = determineReactionOccurrences(firingBounds, violatingTimes, firstViolatingIndex); // Step 4
        
        // if updating the time violates the finite end time value, terminate the simulation
        if (currentTime + violatingTimes[firstViolatingIndex] > endValue && !endInfinity) {
            break;
        }

        updateTime(violatingTimes[firstViolatingIndex]); // Step 5 effecting the leap: time update
        //cout << "TIME ADD: " << violatingTimes[firstViolatingIndex] << endl;
        updateState(reactionOccurrences); // Step 5 effecting the leap: state update

        // record the updated time only if the finalOnly flag is false
        if (!finalOnly) {
            allTimes.push_back(currentTime);
            allStates.push_back(currentState);
        }
        currentIteration++; // update iteration
        props = calculatePropensities(); // Step 1a, sets nonePossible
    }    
}
