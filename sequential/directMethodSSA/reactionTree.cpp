#include "reactionTree.h"

using namespace std;

double reactionTree::calculatePropensity(double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int> > reactants){
    double propensity = reactionRate;
    for (pair<int, int> reactant: reactants) {
        for (int i=0; i<reactant.second; i++) {
            propensity *= (moleculeAmounts[reactant.first] - i);
        }
    }
    return propensity;
}

reactionTree::reactionTree(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector) {
    cout << "Beginning Creation of Reaction Tree" << endl;
    long numReactions = reactantsVector.size();
    reactionTreeArray = new reactionNode[numReactions];
    reactionTreeArray[0].parent = -1;
    for (int i=0; i<numReactions; i++) {
        reactionTreeArray[i].propensity = calculatePropensity(reactionRates[i], moleculeAmounts, reactantsVector[i]);
        reactionTreeArray[i].leftSum = 0;
        reactionTreeArray[i].rightSum = 0;
    }
    for (int i=0; i<numReactions; i++) {
        if (i * 2 + 1 < numReactions) {
            reactionTreeArray[i].leftChild = i * 2 + 1;
            reactionTreeArray[i * 2 + 1].parent = i;
        }
        else {
            reactionTreeArray[i].leftChild = 0;
        }
        if (i * 2 + 2 < numReactions) {
            reactionTreeArray[i].rightChild = i * 2 + 2;
            reactionTreeArray[i * 2 + 2].parent = i;
        }
        else {
            reactionTreeArray[i].rightChild = 0;
        }
    }

    for (long i=numReactions-1; i>0; i--) {
        double subTotalPropensity = reactionTreeArray[i].propensity + reactionTreeArray[i].rightSum + reactionTreeArray[i].leftSum;
        if (i == reactionTreeArray[reactionTreeArray[i].parent].leftChild) {
            reactionTreeArray[reactionTreeArray[i].parent].leftSum += subTotalPropensity;
        }
        else {
            reactionTreeArray[reactionTreeArray[i].parent].rightSum += subTotalPropensity;
        }

    }
}

// Search for node in reaction tree based on uniform RV
int reactionTree::searchForNode(double RV) {
    int currentIndex = 0;
    reactionNode checkNode = reactionTreeArray[0];
    double leftSumTotal = reactionTreeArray[0].leftSum;
    double totalPropensity = reactionTreeArray[0].leftSum + reactionTreeArray[0].rightSum + reactionTreeArray[0].propensity;
    while(RV < (leftSumTotal/totalPropensity) || RV > ((leftSumTotal+checkNode.propensity)/totalPropensity)){
        if (RV < (leftSumTotal/totalPropensity)) {
            currentIndex = currentIndex*2 + 1;
            leftSumTotal -= checkNode.leftSum;
            checkNode = reactionTreeArray[checkNode.leftChild];
            leftSumTotal += checkNode.leftSum;

        } else {
            currentIndex = currentIndex*2 + 2;
            leftSumTotal += checkNode.propensity;
            checkNode = reactionTreeArray[checkNode.rightChild];
            leftSumTotal += checkNode.leftSum;
        }
    }
    return currentIndex;
}

//update its propensity and every subsequent parent's propensities
void reactionTree::updatePropensity(int index, double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int> > reactants){
    int currentIndex = index;
    double newPropensity = calculatePropensity(reactionRate, moleculeAmounts, reactants); //recalculate propensity for that index
    double propensityChange = newPropensity - reactionTreeArray[currentIndex].propensity;
    reactionTreeArray[currentIndex].propensity = newPropensity;
    while (reactionTreeArray[currentIndex].parent != -1) {
        if (currentIndex == reactionTreeArray[reactionTreeArray[currentIndex].parent].leftChild) {
            reactionTreeArray[reactionTreeArray[currentIndex].parent].leftSum += propensityChange;
        }
        else {
            reactionTreeArray[reactionTreeArray[currentIndex].parent].rightSum += propensityChange;
        }
        currentIndex = reactionTreeArray[currentIndex].parent;
    }
}
