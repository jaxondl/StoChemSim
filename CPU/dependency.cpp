#include "dependency.h"


using namespace std;

void dependency::testFunction() {
    cout << "This is a test function for dependency.cpp" << endl;
}

dependency::dependency(string moleculeTypes[], int moleculeAmounts[], vector<vector<pair<int, int>>> stateChangeArray, vector<vector<pair<int,int>>> reactantsArray){
    cout << "Beginning Creation of Dependency Graph" << endl;
    // DEFINITION 1: Reactants(p) and products(p) are reactants and prods of reaction p. e.g. Reactants(1) = {a,b}
    // DEFINITION 2: DependsOn(a-mu), where a-mu is the propensity of chosen reaction, is the set of substances that affect its value. i.e. Reactants(mu)
    // DEFINITION 3: Affects(mu) is set of substances that change in number when reaction mu executes. i.e. Reactants(mu) UNION Products(mu)
    /** DEFINITION 4 (DEPENDENCY GRAPH): Directed graph where vertex set=R (all reactions),
     * directed edge FROM Vi to Vj IFF Affects(Vi) INTERSECTION DependsOn(a-Vj) is NOT an empty set.
     * i.e. at least one of the reactants and products of Vi is shared with the reactants of Vj.
    */

    int numReactions = stateChangeArray.size();
    vector<set<string>> dependsOn(numReactions);
    vector<set<string>> affects(numReactions);

    for(int i = 0; i < numReactions; i++){
        for(pair<int, int> element : reactantsArray[i]){
            dependsOn[i].insert(moleculeTypes[element.first]);
        }
        for(pair<int, int> element : stateChangeArray[i]){
            affects[i].insert(moleculeTypes[element.first]);
        }
    }

    vector<vector<int>> dummyGraph(numReactions);

    //then find intersection
    for(int i = 0; i < numReactions; i++){
        dummyGraph[i].push_back(i);
        for(int j = 0; j < numReactions; j++){
            if(j != i && intersects(affects[i], dependsOn[j])){
                dummyGraph[i].push_back(j);
            }
        }
    }

    dependencyGraph = dummyGraph;
}

bool dependency::intersects(set<string> set1, set<string> set2){
    for(string s : set1){
        cout << s << " ";
    }
    cout << endl;
    for(string s: set2){
        cout << s << " ";
    }
    cout << endl;

    //begin intersection algorithm
    set<string>::iterator iter1 = set1.begin();
    set<string>::iterator iter2 = set2.begin();
    while(iter1 != set1.end() && iter2 != set2.end()){
        if((*iter1).compare(*iter2) < 0){ //that means *iter1 is lexicographically smaller than *iter2
            iter1++;
        }
        else if((*iter1).compare(*iter2) > 0){ //that means *iter2 is lexicographically smaller than *iter1
            iter2++;
        }
        else{
            return true;
        }
    }
    return false;
}