#include "reactionTree.h"

using namespace std;

// Propesnity calculation for each reaction = k * (amount(reactantMolecule1)) * (amount(reactantMolecule1)-1)...(amount(reactantMolecule1)-reactantCoefficient + 1) * ... continue for every reactant
double reactionTree::calculatePropensity(double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int> > reactants){
    double propensity = reactionRate;
    for (pair<int, int> reactant: reactants) { 
        for (int i=0; i<reactant.second; i++) {
            propensity *= (moleculeAmounts[reactant.first] - i);
        }
    }
    return propensity;
}

// Create the reaction tree
reactionTree::reactionTree(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector) {
    cout << "Beginning Creation of Reaction Tree" << endl;
    int numReactions = reactantsVector.size();
    reactionTreeArray = new reactionNode[numReactions];
    reactionTreeArray[0].parent = -1;
    // Calculate propensity for each reaction
    for (int i=0; i<numReactions; i++) {
        reactionTreeArray[i].propensity = calculatePropensity(reactionRates[i], moleculeAmounts, reactantsVector[i]);
        reactionTreeArray[i].leftSum = 0;
        reactionTreeArray[i].rightSum = 0;
    }
    // Assign left and right child indices for each tree node beginning at the roott
    for (int i=0; i<numReactions; i++) {
        if (i * 2 + 1 < numReactions) { 
            reactionTreeArray[i].leftChild = i * 2 + 1;
            reactionTreeArray[i * 2 + 1].parent = i;
        }
        else {
            reactionTreeArray[i].leftChild = -1;
        }
        if (i * 2 + 2 < numReactions) {
            reactionTreeArray[i].rightChild = i * 2 + 2;
            reactionTreeArray[i * 2 + 2].parent = i;
        }
        else {
            reactionTreeArray[i].rightChild = -1;
        }
    }

    // Calculating left and right sums for all tree nodes starting from the leaves (end of array to the beginning)
    for (int i=numReactions-1; i>0; i--) {
        double subTotalPropensity = reactionTreeArray[i].propensity + reactionTreeArray[i].rightSum + reactionTreeArray[i].leftSum; // current node's subtree propensity sum
        // Add subTotalPropensity to parent's leftSum or rightSum depending on if the child is the left or right child of its parent
        if (i == reactionTreeArray[reactionTreeArray[i].parent].leftChild) {
            reactionTreeArray[reactionTreeArray[i].parent].leftSum += subTotalPropensity;
        }
        else {
            reactionTreeArray[reactionTreeArray[i].parent].rightSum += subTotalPropensity;
        }
    }
}

// Search for node in reaction tree based on sampled uniform RV
int reactionTree::searchForNode(double RV) {
    int currentIndex = 0; // index of checkNode in the reactionTree array, start with root node
    reactionNode checkNode = reactionTreeArray[0]; // node object of the reaction to be checked
    double leftSumTotal = reactionTreeArray[0].leftSum;
    double totalPropensity = reactionTreeArray[0].leftSum + reactionTreeArray[0].rightSum + reactionTreeArray[0].propensity;
    while(RV < (leftSumTotal/totalPropensity) || RV > ((leftSumTotal+checkNode.propensity)/totalPropensity)) { // if the sampled RV is not in the current node's propensity range, continue searching
        if (RV < (leftSumTotal/totalPropensity)) { // if the sampled RV is in the left subtree of the current node, update the current node to the left child
            currentIndex = currentIndex*2 + 1;
            leftSumTotal -= checkNode.leftSum;
            if (checkNode.leftChild == -1){
                return -1;
                cout << "PROBLEM" << endl;
                cout << "RV value: " << RV << endl;
                cout << "Left Sum Total" << leftSumTotal << " Propensity: " << checkNode.propensity << endl;
                break;
            }
            checkNode = reactionTreeArray[checkNode.leftChild];
            leftSumTotal += checkNode.leftSum;
        } else { // if the sampled RV is in the right subtree of the current node, update the current node to the right child
            currentIndex = currentIndex*2 + 2;
            leftSumTotal += checkNode.propensity;
            if (checkNode.rightChild == -1){
                return -1;
                cout << "PROBLEM" << endl;
                cout << "RV value: " << RV << endl;
                cout << "Left Sum Total" << leftSumTotal << " Propensity: " << checkNode.propensity << endl;
                break;
            }
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
    double propensityChange = newPropensity - reactionTreeArray[currentIndex].propensity; // determine change in propensity
    reactionTreeArray[currentIndex].propensity = newPropensity;
    while (reactionTreeArray[currentIndex].parent != -1) { // update parents' left or right sums depending on if the child is the left or right child of its parent until you hit the root node
        if (currentIndex == reactionTreeArray[reactionTreeArray[currentIndex].parent].leftChild) {
            reactionTreeArray[reactionTreeArray[currentIndex].parent].leftSum += propensityChange;
        }
        else {
            reactionTreeArray[reactionTreeArray[currentIndex].parent].rightSum += propensityChange;
        }
        currentIndex = reactionTreeArray[currentIndex].parent; // move up the tree to next parent
    }
}
