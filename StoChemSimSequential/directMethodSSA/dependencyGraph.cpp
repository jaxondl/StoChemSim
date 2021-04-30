#include "dependencyGraph.h"

using namespace std;

dependencyGraph::dependencyGraph(vector<vector<pair<int, int> > > stateChangeVector, vector<vector<pair<int,int> > > reactantsVector, vector<int> moleculeAmounts){
    cout << "Beginning Creation of Dependency Graph" << endl;
    // DEFINITION 1: Reactants(p) and products(p) are reactants and prods of reaction p. e.g. Reactants(1) = {a,b}
    // DEFINITION 2: DependsOn(a-mu), where a-mu is the propensity of chosen reaction, is the set of substances that affect its value. i.e. Reactants(mu)
    // DEFINITION 3: Affects(mu) is set of substances that change in number when reaction mu executes. i.e. Reactants(mu) UNION Products(mu)
    /** DEFINITION 4 (DEPENDENCY GRAPH): Directed graph where vertex set=R (all reactions),
     * directed edge FROM Vi to Vj IFF Affects(Vi) INTERSECTION DependsOn(a-Vj) is NOT an empty set.
     * i.e. at least one of the reactants and products of Vi is shared with the reactants of Vj.
    */

    // A+B=C reaction 1
    // dependson(1) = [A, B]
    // affects(1) = [A,B,C]

    // draw an edge from reaction 1 to reaction 2 IFF dependsOn(2) has something in affects(1)

    //no need for an affects vector if just going to use state change molecules
    int numMolecules = moleculeAmounts.size();
    vector<vector<int> > dependsOn(numMolecules);  // each molecule index will have a set of the reaction indices of which have the molecule as a reactant
    for(int i = 0; i < reactantsVector.size(); i++){
        for(pair<int, int> element : reactantsVector[i]){
            dependsOn[element.first].push_back(i); // insert the reaction index in to the molecule's depends on
        }
    }

    vector<vector<int> > dummyGraph(stateChangeVector.size());
    int size = stateChangeVector.size() * stateChangeVector.size();
    int8_t* checkAdded = new int8_t[size];
    memset(checkAdded, 0, size);
    for(int i = 0; i < stateChangeVector.size(); i++){
        int index_base = stateChangeVector.size() * i;
        for(pair<int, int> element : stateChangeVector[i]){
            for(auto dep : dependsOn[element.first]){
                int index = index_base + dep;
                if (checkAdded[index] == 0){
                    checkAdded[index] = 1;
                    dummyGraph[i].push_back(dep);
                }
            }
        }
    }

    dependencyGraphStructure = dummyGraph; // assign the "dummy" graph to the actual dependencyGraphStructure within the class
}

bool dependencyGraph::intersects(set<int> set1, set<int> set2){
    set<int>::iterator iter1 = set1.begin(); // create two iterators for both sets in question
    set<int>::iterator iter2 = set2.begin();
    while(iter1 != set1.end() && iter2 != set2.end()){ // while either set has not been completely traversed
        if(*iter1 < *iter2){ // this means *iter1 is lexicographically smaller than *iter2 for string; in this case that means the index is smaller
            iter1++; // thus, look at the next value in set 1
        }
        else if(*iter1 > *iter2){ // this means *iter2 is lexicographically smaller than *iter1 instead
            iter2++; // instead, look at the next value in set 2
        }
        else{
            return true; // otherwise, the two values must be lexicographically identical. A match, or intersection, is found
        }
    }
    return false; // if either set has already been traversed and no values are identical, then there is no intersection
}

vector<int> dependencyGraph::getDependentReactions(int reactionIndex){
    return dependencyGraphStructure[reactionIndex]; // returns the designated vector within the dependency graph for that reaction index
}
