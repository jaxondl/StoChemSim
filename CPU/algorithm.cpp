#include "algorithm.h"

using namespace std;

algorithm::algorithm(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double t_end){

    reactionTree = new tree(moleculeAmounts, reactionRates, reactantsVector);
    dependencyGraph = new dependency(stateChangeVector, reactantsVector);

    reactionRates = reactionRates; // k reaction constants
    reactantsVector = reactantsVector;
    stateChangeVector = stateChangeVector;
    currentState = moleculeAmounts;
    currentTime = 0;
    t_end = t_end;
}

double algorithm::getUniformRandomVariable(){
    random_device rd; // random seed
    mt19937 gen(rd()); // mersenne twister engine
    uniform_real_distribution<> dis(0.0, 1.0);
    double RV = dis(gen);
    return RV;
}

void algorithm::updateTime(double timeUntilNextReaction){
    currentTime += timeUntilNextReaction;
}

void algorithm::updateState(vector<vector<pair<int, int>>> stateChangeVector, int reactionIndex) {
    vector<pair<int, int>> chosenReactionChange = stateChangeVector[reactionIndex];
    for(pair<int, int> p: chosenReactionChange){
        currentState[p.first] += p.second; //update the state change at the associated index
    }
}

double algorithm::getTimeUntilNextReaction(double propensity) {
    random_device rd;
    mt19937 gen(rd());
    exponential_distribution<> dis(propensity); //exponential approximation of propensity
    double RV = dis(gen);
    return RV;
}

double algorithm::getTotalPropensity(){
    tree::ReactionNode root = reactionTree->ReactionTreeArray[0];
    double totalPropensity = root.propensity + root.leftSum + root.rightSum;
    return totalPropensity;
}

void algorithm::start(){
    /**Step 1:
     */
    while (currentTime < t_end && getTotalPropensity() != 0){
        /**Step 1: calculate time until next reaction by doing getTimeUntilNectreaction(getTotalPropensity)
         * step 2: update time += time until next reaction
         * step 3: get uniform
         * Step 4: refer to tree and pick the bucket (reaction)
         * step 5: send reaction number to depedency graph and get the dependent reactions
         * step 6: for loop and update props in tree w fxn
         * step 7: update state using state change vector
         * (repeat) 
        */
        double timeUntilNextReaction = getTimeUntilNextReaction(getTotalPropensity());
        updateTime(timeUntilNextReaction);
        double uniformRV = getUniformRandomVariable();
        int reactionIndex = reactionTree->searchForNode(uniformRV);
        vector<int> dependentReactionIndices = dependencyGraph->getDependentReactions(reactionIndex);
        for(int reaction: dependentReactionIndices){
            reactionTree->updatePropensity(reaction, reactionRates[reaction], currentState, reactantsVector[reaction]);
        }
        updateState(stateChangeVector, reactionIndex);
    }
}