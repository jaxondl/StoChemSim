#include "algorithm.h"

using namespace std;

algorithm::algorithm(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double tEnd){

    reactionTree = new tree(moleculeAmounts, reactionRates, reactantsVector);
    dependencyGraph = new dependency(stateChangeVector, reactantsVector);

    this->reactionRates = reactionRates; // k reaction constants
    this->reactantsVector = reactantsVector;
    this->stateChangeVector = stateChangeVector;
    this->currentState = moleculeAmounts;
    this->currentTime = 0;
    this->t_end = tEnd;
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
        cout << "Current amount: " << currentState[p.first] << " state change: " << p.second << endl;
        currentState[p.first] += p.second; //update the state change at the associated index
    }

}

double algorithm::getTimeUntilNextReaction(double propensity) {
    random_device rd;
    mt19937 gen(rd());
    cout << "HELLO THERE" << endl;
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
    while (getTotalPropensity() > 0.001){
        double timeUntilNextReaction = 0.37; // getTimeUntilNextReaction(getTotalPropensity());
        updateTime(timeUntilNextReaction);
        double uniformRV = getUniformRandomVariable();
        int reactionIndex = reactionTree->searchForNode(uniformRV);
        updateState(stateChangeVector, reactionIndex);
        vector<int> dependentReactionIndices = dependencyGraph->getDependentReactions(reactionIndex);
        for(int reaction: dependentReactionIndices){
            reactionTree->updatePropensity(reaction, reactionRates[reaction], currentState, reactantsVector[reaction]);
        }
        for (int molecule: currentState){
            cout << molecule << " ";
        }
        cout << endl;
        cout << getTotalPropensity() << endl;
    }
    double d = 0.0;
    for (int i = 0; i< reactionRates.size(); i++){
        cout << "Propensity at index: " << i << " " << this->reactionTree->ReactionTreeArray[i].propensity << endl;
        d += this->reactionTree->ReactionTreeArray[i].propensity;
    }
    cout << d;
}
/**Step 1: calculate time until next reaction by doing getTimeUntilNectreaction(getTotalPropensity)
         * step 2: update time += time until next reaction
         * step 3: get uniform
         * Step 4: refer to tree and pick the bucket (reaction)
         * step 5: send reaction number to depedency graph and get the dependent reactions
         * step 6: for loop and update props in tree w fxn
         * step 7: update state using state change vector
         * (repeat) 
        */