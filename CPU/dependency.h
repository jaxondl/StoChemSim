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
    int testInt;
public:
    vector<vector<int>> dependencyGraph;

    dependency(string moleculeTypes[], int moleculeAmounts[], vector<vector<pair<int, int>>> stateChangeArray, vector<vector<pair<int,int>>> reactantsArray);

    bool intersects(set<string> set1, set<string> set2);

    void testFunction();
};


#endif //DEPENDENCY_H
