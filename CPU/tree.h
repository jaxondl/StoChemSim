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
    double calculatePropensity(double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int>> reactants);
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
    tree(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector);
    int searchForNode(double RV);
    void updatePropensity(int index, double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int>> reactants);
};


#endif //TREE_H
