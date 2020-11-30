#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <string>
#include <iterator>
#include <set>
#include <algorithm>
#include <vector>
#include <utility>

using namespace std;

class tree {
private:
    double calculatePropensity(double reactionRate, int moleculeAmounts[], vector<pair<int, int>> reactants);
public:
    struct ReactionNode {
            double propensity;
            double leftSum;
            double rightSum;
            int leftChild;
            int rightChild;
            int parent;
    };
    ReactionNode* ReactionTreeArray;
    tree(int numReactions, int moleculeAmounts[], double reactionRates[], vector<pair<int, int>> ReactantsArray[]);
    int searchForNode(double RV);
    void updatePropensity(int index, double reactionRate, int moleculeAmounts[], vector<pair<int, int>> reactants);
};


#endif //TREE_H
