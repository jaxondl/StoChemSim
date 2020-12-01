#ifndef DEPENDENCY_H
#define DEPENDENCY_H

#include <iostream>
#include <string>
#include <iterator>
#include <set>
#include <algorithm>
#include <vector>
#include <utility>

using namespace std;

class dependency {
private:
    vector<vector<int>> dependencyGraph; 
public:
    dependency(vector<vector<pair<int, int>>> stateChangeVector, vector<vector<pair<int,int>>> reactantsVector);

    bool intersects(set<int> set1, set<int> set2);

    vector<int> getDependentReactions(int reactionIndex);
};


#endif //DEPENDENCY_H
