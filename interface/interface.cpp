/* Include required headers */
#include <cstdint>
#include <vector>
#include "WolframLibrary.h"
#include "WolframNumericArrayLibrary.h"

#include <iostream>
#include <string>
#include <iterator>
#include <set>
#include <random>
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

dependencyGraph::dependencyGraph(vector<vector<pair<int, int>>> stateChangeVector, vector<vector<pair<int,int>>> reactantsVector){
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

// ========== Dependency Graph ===========

class reactionTree {
private:
    double calculatePropensity(double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int>> reactants);
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
    reactionTree(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector);
    int searchForNode(double RV);
    void updatePropensity(int index, double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int>> reactants);
};

double reactionTree::calculatePropensity(double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int>> reactants){
    double propensity = reactionRate;
    for (pair<int, int> reactant: reactants) {
        for (int i=0; i<reactant.second; i++) {
            propensity *= (moleculeAmounts[reactant.first] - i);
        }
    }
    if(reactants.size() == 0){
        propensity = 0;
    }
    return propensity;
}

reactionTree::reactionTree(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector) {
    long numReactions = reactantsVector.size();
    reactionTreeArray = new reactionNode[numReactions];
    reactionTreeArray[0].parent = -1;
    for (int i=0; i<numReactions; i++) {
        reactionTreeArray[i].propensity = calculatePropensity(reactionRates[i], moleculeAmounts, reactantsVector[i]);
        reactionTreeArray[i].leftSum = 0;
        reactionTreeArray[i].rightSum = 0;
    }
    for (int i=0; i<numReactions; i++) {
        if (i * 2 + 1 < numReactions) {
            reactionTreeArray[i].leftChild = i * 2 + 1;
            reactionTreeArray[i * 2 + 1].parent = i;
        }
        else {
            reactionTreeArray[i].leftChild = 0;
        }
        if (i * 2 + 2 < numReactions) {
            reactionTreeArray[i].rightChild = i * 2 + 2;
            reactionTreeArray[i * 2 + 2].parent = i;
        }
        else {
            reactionTreeArray[i].rightChild = 0;
        }
    }

    for (long i=numReactions-1; i>0; i--) {
        double subTotalPropensity = reactionTreeArray[i].propensity + reactionTreeArray[i].rightSum + reactionTreeArray[i].leftSum;
        if (i == reactionTreeArray[reactionTreeArray[i].parent].leftChild) {
            reactionTreeArray[reactionTreeArray[i].parent].leftSum += subTotalPropensity;
        }
        else {
            reactionTreeArray[reactionTreeArray[i].parent].rightSum += subTotalPropensity;
        }

    }
}

// Search for node in reaction tree based on uniform RV
int reactionTree::searchForNode(double RV) {
    int currentIndex = 0;
    reactionNode checkNode = reactionTreeArray[0];
    double leftSumTotal = reactionTreeArray[0].leftSum;
    double totalPropensity = reactionTreeArray[0].leftSum + reactionTreeArray[0].rightSum + reactionTreeArray[0].propensity;
    while(RV < (leftSumTotal/totalPropensity) || RV > ((leftSumTotal+checkNode.propensity)/totalPropensity)){
        if (RV < (leftSumTotal/totalPropensity)) {
            currentIndex = currentIndex*2 + 1;
            leftSumTotal -= checkNode.leftSum;
            checkNode = reactionTreeArray[checkNode.leftChild];
            leftSumTotal += checkNode.leftSum;

        } else {
            currentIndex = currentIndex*2 + 2;
            leftSumTotal += checkNode.propensity;
            checkNode = reactionTreeArray[checkNode.rightChild];
            leftSumTotal += checkNode.leftSum;
        }
    }
    return currentIndex;
}

//update its propensity and every subsequent parent's propensities
void reactionTree::updatePropensity(int index, double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int>> reactants){
    int currentIndex = index;
    double newPropensity = calculatePropensity(reactionRate, moleculeAmounts, reactants); //recalculate propensity for that index
    double propensityChange = newPropensity - reactionTreeArray[currentIndex].propensity;
    reactionTreeArray[currentIndex].propensity = newPropensity;
    while (reactionTreeArray[currentIndex].parent != -1) {
        if (currentIndex == reactionTreeArray[reactionTreeArray[currentIndex].parent].leftChild) {
            reactionTreeArray[reactionTreeArray[currentIndex].parent].leftSum += propensityChange;
        }
        else {
            reactionTreeArray[reactionTreeArray[currentIndex].parent].rightSum += propensityChange;
        }
        currentIndex = reactionTreeArray[currentIndex].parent;
    }
}

// ========== Reaction Tree ===========

class directMethodSSA {
private:
    reactionTree* reactionTree; 
    dependencyGraph* dependencyGraph;

    vector<double> reactionRates;
    vector<vector<pair<int, int>>> reactantsVector;
    vector<vector<pair<int, int>>> stateChangeVector;
    vector<int> currentState;
    vector<vector<int>> allStates;
    vector<double> allTimes;
    double currentTime;
    double tEnd;

    double getUniformRandomVariable();
    double getTimeUntilNextReaction(double propensity);
    void updateTime(double timeUntilNextReaction);
    void updateState(vector<vector<pair<int, int>>> stateChangeVector, int reactionIndex);
    double getTotalPropensity();

public:
    directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double t_end);
    vector<vector<int>> getAllStates();
    vector<double> getAllTimes();
    void start();
};

directMethodSSA::directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int>>> reactantsVector, vector<vector<pair<int, int>>> stateChangeVector, double tEnd){
    this->reactionTree = new class reactionTree(moleculeAmounts, reactionRates, reactantsVector);
    this->dependencyGraph = new class dependencyGraph(stateChangeVector, reactantsVector);
    this->allStates.push_back(moleculeAmounts);
    this->allTimes.push_back(0);
    this->reactionRates = reactionRates; // k reaction constants
    this->reactantsVector = reactantsVector;
    this->stateChangeVector = stateChangeVector;
    this->currentState = moleculeAmounts;
    this->currentTime = 0;
    this->tEnd = tEnd;
}

double directMethodSSA::getUniformRandomVariable(){
    random_device rd; // random seed
    mt19937 gen(rd()); // mersenne twister engine
    uniform_real_distribution<> dis(0.0, 1.0);
    double RV = dis(gen);
    return RV;
}

void directMethodSSA::updateTime(double timeUntilNextReaction){
    currentTime += timeUntilNextReaction;
}

void directMethodSSA::updateState(vector<vector<pair<int, int>>> stateChangeVector, int reactionIndex) {
    vector<pair<int, int>> chosenReactionChange = stateChangeVector[reactionIndex];
    for(pair<int, int> p: chosenReactionChange){
        currentState[p.first] += p.second; //update the state change at the associated index
    }

}

double directMethodSSA::getTimeUntilNextReaction(double propensity) {
    random_device rd;
    mt19937 gen(rd());
    exponential_distribution<> dis(propensity); //exponential approximation of propensity
    double RV = dis(gen);
    return RV;
}

double directMethodSSA::getTotalPropensity(){
    reactionTree::reactionNode root = reactionTree->reactionTreeArray[0];
    double totalPropensity = root.propensity + root.leftSum + root.rightSum;
    return totalPropensity;
}

vector<vector<int>> directMethodSSA::getAllStates(){
    return allStates;
}

vector<double> directMethodSSA::getAllTimes(){
    return allTimes;
}

void directMethodSSA::start(){
    while (getTotalPropensity() > 0.001 && currentTime < tEnd){
        double timeUntilNextReaction = getTimeUntilNextReaction(getTotalPropensity());
        updateTime(timeUntilNextReaction);
        allTimes.push_back(currentTime);
        double uniformRV = getUniformRandomVariable();
        int reactionIndex = reactionTree->searchForNode(uniformRV);
        updateState(stateChangeVector, reactionIndex);
        allStates.push_back(currentState);
        vector<int> dependentReactionIndices = dependencyGraph->getDependentReactions(reactionIndex);
        for(int reaction: dependentReactionIndices){
            reactionTree->updatePropensity(reaction, reactionRates[reaction], currentState, reactantsVector[reaction]);
        }
    }
    // If we pass tEnd, remove the last time and state before sending 
    if (currentTime > tEnd) {
        allTimes.pop_back();
        allStates.pop_back();
    }
}

// ========== Direct SSA ===========

// ========= backend code end =======

/* Return the version of Library Link */
DLLEXPORT mint WolframLibrary_getVersion() { return WolframLibraryVersion; }

/* Initialize Library */
DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData) {
	return LIBRARY_NO_ERROR;
}

/* Uninitialize Library */
DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData) {
	return;
}

/* convert a 2-dimentional NumericArray to a 2-dimentional vector */
template <typename Tin, typename Tout>
static vector<vector<int64_t>> numericMatrixtoVector(const void *in0, mint const *dims) {
	const Tin *in = static_cast<const Tin *>(in0);
	vector<vector<Tout>> out;
	mint row = dims[0];
	mint col = dims[1];
	for (mint i = 0; i < row; i++) {
		vector<T> out_row;
		for (mint j = 0; j < col; j++) {
			out_row.push_back(in[i*col + j]);
		}
		out.push_back(out_row);
	}
	return out;
}

/* convert a 1-dimentional NumericArray to a 1-dimentional vector */
template <typename Tin, typename Tout>
static vector<Tout> numericArraytoVector(const void *in0, mint const length) {
	const Tin *in = static_cast<const Tin *>(in0);
    vector<Tout> out;
    for (mint i = 0; i < length; i++) {
        out.push_back(in[i]);
    }
	return out;
}

/* convert a 1-dimentional vector to a 1-dimentional NumericArray */
template <typename T>
static void vectortoNumericArray(void *Mout0, vector<T> out) {
	T *Mout = static_cast<T *>(Mout0);
	for (int64_t i = 0; i < out.size(); i++) {
		Mout[i] = out[i];
	}
}

template <typename Tin, Tout>
static void matrixtoNumericArray(void *Mout0, vector<<vector<Tin>>> out) {
	Tout *Mout = static_cast<Tout *>(Mout0);
	int64_t row = out.size();
	int64_t col = out[0].size();
	for (int64_t i = 0; i < row; i++) {
		for (int64_t j = 0; j < col; j++) {
			Mout[i*col + j] = (Tout)out[i][j];
		}
	}
}

template <typename T1, typename T2>
static void reactantsAndStateChangeArrayConstruction(const int64_t *reactIn, const int64_t *prodIn, vector<vector<pair<T1, T2>>>& reactantsArray, vector<vector<pair<T1, T2>>>& stateChangeArray) {
	for (mint i = 0; i < reactionCount; i++) {
		vector<pair<T1, T2>> reactantsArray_row;
		vector<pair<T1, T2>> stateChangeArray_row;
		for (mint j = 0; j < moleculeCount; j++) {
			T1 index = (T1)j;
			T2 in = (T2)reactIn[i*moleculeCount + j];
			T2 out = (T2)prodIn[i*moleculeCount + j];
			if(in > 0) {
				reactantsArray_row.push_back(pair<T1, T2>(index, in));
				stateChangeArray_row.push_back(pair<T1, T2>(index, -in));
			}
			if(out > 0) {
				stateChangeArray_row.push_back(pair<T1, T2>(index, out));
			}
		}
        reactantsArray.push_back(reactantsArray_row);
		stateChangeArray.push_back(stateChangeArray_row);
    }
}

// ****** gloabl storage for returns *******
vector<vector<int>> allStates;
vector<double> allTimes;

/* CRN SSA main function */
EXTERN_C DLLEXPORT int CRN_SSA(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument res) {
	// debug setup
	int err = LIBRARY_FUNCTION_ERROR;
	WolframNumericArrayLibrary_Functions naFuns = libData->numericarrayLibraryFunctions;

	// reused local varibles setup
	void *data_out = NULL;

	// convert initCounts
	MNumericArray MinitCounts = MArgument_getMNumericArray(Args[0]);
	void *data_in = naFuns->MNumericArray_getData(MinitCounts);
	mint length = naFuns->MNumericArray_getFlattenedLength(MinitCounts);
	vector<int> moleculeAmounts = numericArraytoVector<int, int>(data_in, length);

	// convert reactantsArray & stateChangeArray
	MNumericArray MreactCounts = MArgument_getMNumericArray(Args[1]);
	MNumericArray MprodCounts = MArgument_getMNumericArray(Args[2]);
	mint const * dims = naFuns->MNumericArray_getDimensions(MreactCounts);
	mint reactionCount = dims[0];
	mint moleculeCount = dims[1];
	void* MreactCounts_in = naFuns->MNumericArray_getData(MreactCounts);
	void* MprodCounts_in = naFuns->MNumericArray_getData(MprodCounts);

	const int64_t *reactIn = static_cast<const int64_t *>(MreactCounts_in);
	const int64_t *prodIn = static_cast<const int64_t *>(MprodCounts_in);

    vector<vector<pair<int, int>>> reactantsArray;
	vector<vector<pair<int, int>>> stateChangeArray;
    reactantsAndStateChangeArrayConstruction<int, int>(reactIn, prodIn, reactantsArray, stateChangeArray);

	// convert rates
	MNumericArray Mrates = MArgument_getMNumericArray(Args[3]);
	data_in = naFuns->MNumericArray_getData(Mrates);
	length = naFuns->MNumericArray_getFlattenedLength(Mrates);
	vector<double> kValues = numericArraytoVector<double, double>(data_in, length);

	// convert tEnd
	mreal tEnd = MArgument_getReal(Args[4]);

	// CRN SSA process: pass everything to backend
	process = new directMethodSSA(
					moleculeAmounts,
					kValues,
					reactantsArray,
					stateChangeArray,
					(double)tEnd);
	process::start();

	allStates = process::getAllStates();
	allTimes = process::getAllTimes();
	
	return LIBRARY_NO_ERROR;
}


EXTERN_C DLLEXPORT int getTimes(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument res) {
	// debug setup
	int err = LIBRARY_FUNCTION_ERROR;
	WolframNumericArrayLibrary_Functions naFuns = libData->numericarrayLibraryFunctions;

	// reused local varibles setup
	void *data_in = NULL, *data_out = NULL;
	mint length;
	mint const *dims;

	// output setup
	MNumericArray Mout;
	int64_t out_size = allTimes.size();
	const mint *dims_out = &out_size;
	err = naFuns->MNumericArray_new(MNumericArray_Type_Real64, 1, dims_out, &Mout);
	if (err != 0) {
		goto cleanup;
	}
	data_out = naFuns->MNumericArray_getData(Mout);
	if (data_out == NULL) {
		goto cleanup;
	}
	
	// convert output to a NumericArray
	vectortoNumericArray<mreal>(data_out, allTimes);

	// pass the result back
	MArgument_setMNumericArray(res, Mout);
	return LIBRARY_NO_ERROR;

	cleanup:
	naFuns->MNumericArray_free(Mout);
	return err;
}

EXTERN_C DLLEXPORT int getStates(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument res) {
	// debug setup
	int err = LIBRARY_FUNCTION_ERROR;
	WolframNumericArrayLibrary_Functions naFuns = libData->numericarrayLibraryFunctions;

	// reused local varibles setup
	void *data_in = NULL, *data_out = NULL;
	mint length;
	mint const *dims;

	// output setup
	MNumericArray Mout;
	int64_t row = allStates.size();
	int64_t col = allStates[0].size();
	int64_t out_size[2] = {row, col};
	const mint *dims_out = &out_size;
	err = naFuns->MNumericArray_new(MNumericArray_Type_Bit64, 2, dims_out, &Mout);
	if (err != 0) {
		goto cleanup;
	}
	data_out = naFuns->MNumericArray_getData(Mout);
	if (data_out == NULL) {
		goto cleanup;
	}
	
	// convert output to a NumericArray
	matrixtoNumericArray<int, mint>(data_out, allStates);

	// pass the result back
	MArgument_setMNumericArray(res, Mout);
	return LIBRARY_NO_ERROR;

	cleanup:
	naFuns->MNumericArray_free(Mout);
	return err;
}
