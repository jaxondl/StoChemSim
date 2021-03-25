#ifndef DEPENDENCYGRAPH_H
#define DEPENDENCYGRAPH_H

#include <iostream>
#include <string>
#include <iterator>
#include <set>
#include <algorithm>
#include <vector>
#include <utility>

using namespace std;

class dependencyGraph {
private:
    vector<vector<int> > dependencyGraphStructure; // the graph structure itself
public:
    dependencyGraph(vector<vector<pair<int, int> > > stateChangeVector, vector<vector<pair<int,int> > > reactantsVector, vector<int> moleculeAmounts); // constructor

    bool intersects(set<int> set1, set<int> set2); // user made function to check for an intersection of 1 or more elements between two sets

    vector<int> getDependentReactions(int reactionIndex); // obtain the dependent reactions for a given reaction via the dependencyGraphStructure
};

#endif
