#include "dependencyGraph.h"

using namespace std;

dependencyGraph::dependencyGraph(vector<vector<pair<int, int>>> stateChangeVector, vector<vector<pair<int,int>>> reactantsVector){
    cout << "Beginning Creation of Dependency Graph" << endl;
    // DEFINITION 1: Reactants(p) and products(p) are reactants and prods of reaction p. e.g. Reactants(1) = {a,b}
    // DEFINITION 2: DependsOn(a-mu), where a-mu is the propensity of chosen reaction, is the set of substances that affect its value. i.e. Reactants(mu)
    // DEFINITION 3: Affects(mu) is set of substances that change in number when reaction mu executes. i.e. Reactants(mu) UNION Products(mu)
    /** DEFINITION 4 (DEPENDENCY GRAPH): Directed graph where vertex set=R (all reactions),
     * directed edge FROM Vi to Vj IFF Affects(Vi) INTERSECTION DependsOn(a-Vj) is NOT an empty set.
     * i.e. at least one of the reactants and products of Vi is shared with the reactants of Vj.
    */

    int numReactions = stateChangeVector.size();
    vector<set<int>> dependsOn(numReactions);
    vector<set<int>> affects(numReactions);

    for(int i = 0; i < numReactions; i++){
        for(pair<int, int> element : reactantsVector[i]){
            dependsOn[i].insert(element.first);
        }
        for(pair<int, int> element : stateChangeVector[i]){
            affects[i].insert(element.first);
        }
    }

    vector<vector<int>> dummyGraph(numReactions);

    for(int i = 0; i < numReactions; i++){ //then find intersection
        dummyGraph[i].push_back(i); 
        for(int j = 0; j < numReactions; j++){
            if(j != i && intersects(affects[i], dependsOn[j])){
                dummyGraph[i].push_back(j);
            }
        }
    }

    dependencyGraphStructure = dummyGraph;
}

bool dependencyGraph::intersects(set<int> set1, set<int> set2){
    set<int>::iterator iter1 = set1.begin();
    set<int>::iterator iter2 = set2.begin();
    while(iter1 != set1.end() && iter2 != set2.end()){
        if(*iter1 < *iter2){ //that means *iter1 is lexicographically smaller than *iter2 for string, but for this case index is smaller
            iter1++;
        }
        else if(*iter1 > *iter2){ //that means *iter2 is lexicographically smaller than *iter1
            iter2++;
        }
        else{
            return true;
        }
    }
    return false;
}

vector<int> dependencyGraph::getDependentReactions(int reactionIndex){
    return dependencyGraphStructure[reactionIndex]; //returns itself as well
}