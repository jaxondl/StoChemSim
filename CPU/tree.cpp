#include "tree.h"

using namespace std;

double tree::calculatePropensity(double reactionRate, int moleculeAmounts[], vector<pair<int, int>> reactants){
    double propensity = reactionRate;
    for (pair<int, int> reactant: reactants) {
        propensity *= pow(moleculeAmounts[reactant.first], reactant.second);
    }
    return propensity;
}

tree::tree(int numReactions, int moleculeAmounts[], double reactionRates[], vector<pair<int, int>> ReactantsArray[]) {
    ReactionTreeArray = new ReactionNode[numReactions];
    //reactionCount = numReactions;
    ReactionTreeArray[0].parent = -1;
    for (int i=0; i<numReactions; i++) {
        ReactionTreeArray[i].propensity = calculatePropensity(reactionRates[i], moleculeAmounts, ReactantsArray[i]);
        ReactionTreeArray[i].leftSum = 0;
        ReactionTreeArray[i].rightSum = 0;
    }
    for (int i=0; i<numReactions; i++) {
        if (i * 2 + 1 < numReactions) {
            ReactionTreeArray[i].leftChild = i * 2 + 1;
            ReactionTreeArray[i * 2 + 1].parent = i;
        }
        else {
            ReactionTreeArray[i].leftChild = 0;
        }
        if (i * 2 + 2 < numReactions) {
            ReactionTreeArray[i].rightChild = i * 2 + 2;
            ReactionTreeArray[i * 2 + 2].parent = i;
        }
        else {
            ReactionTreeArray[i].rightChild = 0;
        }
    }

    for (int i=numReactions-1; i>0; i--) {
        double subTotalPropensity = ReactionTreeArray[i].propensity + ReactionTreeArray[i].rightSum + ReactionTreeArray[i].leftSum;
        if (i == ReactionTreeArray[ReactionTreeArray[i].parent].leftChild) {
            ReactionTreeArray[ReactionTreeArray[i].parent].leftSum += subTotalPropensity;
        }
        else {
            ReactionTreeArray[ReactionTreeArray[i].parent].rightSum += subTotalPropensity;
        }

    }
}

int tree::searchForNode(double RV) {
    int currentIndex = 0;
    ReactionNode checkNode = ReactionTreeArray[0];
    double leftSumTotal = ReactionTreeArray[0].leftSum;
    double totalPropensity = ReactionTreeArray[0].leftSum + ReactionTreeArray[0].rightSum + ReactionTreeArray[0].propensity;
    while(RV < (leftSumTotal/totalPropensity) || RV > ((leftSumTotal+checkNode.propensity)/totalPropensity)){
        if (RV < (leftSumTotal/totalPropensity)) {
            currentIndex = currentIndex*2 + 1;
            leftSumTotal -= checkNode.leftSum;
            checkNode = ReactionTreeArray[checkNode.leftChild];
            leftSumTotal += checkNode.leftSum;

        } else {
            currentIndex = currentIndex*2 + 2;
            leftSumTotal += checkNode.propensity;
            checkNode = ReactionTreeArray[checkNode.rightChild];
            leftSumTotal += checkNode.leftSum;
        }
    }
    return currentIndex;
}

void tree::updatePropensity(int index, double reactionRate, int moleculeAmounts[], vector<pair<int, int>> reactants){ //update its and every subsequent parent's propensities
    int currentIndex = index;
    double newPropensity = calculatePropensity(reactionRate, moleculeAmounts, reactants); //recalculate propensity for that index
    double propensityChange = newPropensity - ReactionTreeArray[currentIndex].propensity;
    ReactionTreeArray[currentIndex].propensity = newPropensity;
    while (ReactionTreeArray[currentIndex].parent != -1) {
        if (currentIndex == ReactionTreeArray[ReactionTreeArray[currentIndex].parent].leftChild) {
            ReactionTreeArray[ReactionTreeArray[currentIndex].parent].leftSum += propensityChange;
        }
        else {
            ReactionTreeArray[ReactionTreeArray[currentIndex].parent].rightSum += propensityChange;
        }
        currentIndex = ReactionTreeArray[currentIndex].parent;
    }
}
