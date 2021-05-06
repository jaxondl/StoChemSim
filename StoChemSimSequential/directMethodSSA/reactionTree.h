#ifndef REACTIONTREE_H
#define REACTIONTREE_H

#include <iostream>
#include <string>
#include <iterator>
#include <set>
#include <algorithm>
#include <vector>
#include <utility>

using namespace std;

class reactionTree {
private:
    double calculatePropensity(double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int> > reactants);
public:
    struct reactionNode {
            double propensity;
            double leftSum;
            double rightSum;
            int leftChild;
            int rightChild;
            int parent;
    };
    reactionNode* reactionTreeArray;
    reactionTree(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector);
    int searchForNode(double RV);
    void updatePropensity(int index, double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int> > reactants);
};

#endif
