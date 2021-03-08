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
    vector<vector<int>> dependencyGraphStructure;
public:
    dependencyGraph(vector<vector<pair<int, int>>> stateChangeVector, vector<vector<pair<int,int>>> reactantsVector);

    bool intersects(set<int> set1, set<int> set2);

    vector<int> getDependentReactions(int reactionIndex);
};

#endif
