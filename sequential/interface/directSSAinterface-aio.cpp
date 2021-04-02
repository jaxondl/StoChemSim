/* Include required headers */
#include <cstdint>
#include <vector>
#include "WolframLibrary.h"
#include "WolframNumericArrayLibrary.h"

/* backend start */
#include <iostream>
#include <string>
#include <iterator>
#include <set>
#include <algorithm>
#include <vector>
#include <utility>
#include <random>

using namespace std;

class dependencyGraph {
private:
    vector<vector<int> > dependencyGraphStructure;
public:
    dependencyGraph(vector<vector<pair<int, int> > > stateChangeVector, vector<vector<pair<int,int> > > reactantsVector);

    bool intersects(set<int> set1, set<int> set2);

    vector<int> getDependentReactions(int reactionIndex);
};

dependencyGraph::dependencyGraph(vector<vector<pair<int, int> > > stateChangeVector, vector<vector<pair<int,int> > > reactantsVector){
    cout << "Beginning Creation of Dependency Graph" << endl;
    // DEFINITION 1: Reactants(p) and products(p) are reactants and prods of reaction p. e.g. Reactants(1) = {a,b}
    // DEFINITION 2: DependsOn(a-mu), where a-mu is the propensity of chosen reaction, is the set of substances that affect its value. i.e. Reactants(mu)
    // DEFINITION 3: Affects(mu) is set of substances that change in number when reaction mu executes. i.e. Reactants(mu) UNION Products(mu)
    /** DEFINITION 4 (DEPENDENCY GRAPH): Directed graph where vertex set=R (all reactions),
     * directed edge FROM Vi to Vj IFF Affects(Vi) INTERSECTION DependsOn(a-Vj) is NOT an empty set.
     * i.e. at least one of the reactants and products of Vi is shared with the reactants of Vj.
    */

    //no need for an affects vector if just going to use state change molecules
    int numMolecules = moleculeAmounts.size();
    vector<set<int>> dependsOn(numMolecules);  // each molecule index will have a set of the reaction indices of which have the molecule as a reactant

    for(int i = 0; i < reactantsVector.size(); i++){
        for(pair<int, int> element : reactantsVector[i]){
            dependsOn[element.first].insert(i); // insert the reaction index in to the molecule's depends on
        }
    }

    vector<vector<int> > dummyGraph(stateChangeVector.size());

    for(int i = 0; i < stateChangeVector.size(); i++){
        for(pair<int, int> element : stateChangeVector[i]){
            for(auto dep : dependsOn[element.first]){
                dummyGraph[i].push_back(dep);
            }
        }
    }

    dependencyGraphStructure = dummyGraph; // assign the "dummy" graph to the actual dependencyGraphStructure within the class
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

double reactionTree::calculatePropensity(double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int> > reactants){
    double propensity = reactionRate;
    for (pair<int, int> reactant: reactants) {
        for (int i=0; i<reactant.second; i++) {
            propensity *= (moleculeAmounts[reactant.first] - i);
        }
    }
    return propensity;
}

reactionTree::reactionTree(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector) {
    cout << "Beginning Creation of Reaction Tree" << endl;
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
void reactionTree::updatePropensity(int index, double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int> > reactants){
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

class directMethodSSA {
private:
    reactionTree* reaction_tree; 
    dependencyGraph* dependency_graph;

    vector<double> reactionRates;
    vector<vector<pair<int, int> > > reactantsVector;
    vector<vector<pair<int, int> > > stateChangeVector;
    vector<int> currentState;
    vector<vector<int> > allStates;
    vector<double> allTimes;
    double currentTime;
    int currentIteration;
    double endValue;
    bool endInfinity;
    bool statesOnly;
    bool finalOnly;
    bool endByIteration;

    double getUniformRandomVariable();
    double getTimeUntilNextReaction(double propensity);
    void updateTime(double timeUntilNextReaction);
    void updateState(vector<vector<pair<int, int> > > stateChangeVector, int reactionIndex);
    double getTotalPropensity();

public:
    directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool statesOnly, bool finalOnly, bool endInfinity, bool endByIteration);
    vector<vector<int> > getAllStates();
    vector<double> getAllTimes();
    vector<int> getCurrentState();
    double getCurrentTime();
    int getCurrentIteration();
    void start();
};

directMethodSSA::directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool statesOnly, bool finalOnly, bool endInfinity, bool endByIteration){
    this->reaction_tree = new class reactionTree(moleculeAmounts, reactionRates, reactantsVector);
    this->dependency_graph = new class dependencyGraph(stateChangeVector, reactantsVector);
    this->allStates.push_back(moleculeAmounts);
    this->allTimes.push_back(0);
    this->reactionRates = reactionRates; // k reaction constants
    this->reactantsVector = reactantsVector;
    this->stateChangeVector = stateChangeVector;
    this->currentState = moleculeAmounts;
    this->currentTime = 0;
    this->currentIteration = 0;
    this->endValue = endValue;
    this->statesOnly = statesOnly;
    this->finalOnly = finalOnly;
    this->endInfinity = endInfinity;
    this->endByIteration = endByIteration;
    if (this->statesOnly && !this->endByIteration) {
        this->endInfinity = true;
    }
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

void directMethodSSA::updateState(vector<vector<pair<int, int> > > stateChangeVector, int reactionIndex) {
    vector<pair<int, int> > chosenReactionChange = stateChangeVector[reactionIndex];
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
    reactionTree::reactionNode root = reaction_tree->reactionTreeArray[0];
    double totalPropensity = root.propensity + root.leftSum + root.rightSum;
    return totalPropensity;
}

vector<vector<int> > directMethodSSA::getAllStates(){
    return allStates;
}

vector<double> directMethodSSA::getAllTimes(){
    return allTimes;
}

vector<int> directMethodSSA::getCurrentState(){
    return currentState;
}

double directMethodSSA::getCurrentTime(){
    return currentTime;
}

int directMethodSSA::getCurrentIteration(){
    return currentIteration;
}

void directMethodSSA::start(){
    while (getTotalPropensity() > 0.001 && ((!endByIteration && (currentTime < endValue || endInfinity)) || (endByIteration && (currentIteration < endValue || endInfinity)))){
        if (!statesOnly) {
            double timeUntilNextReaction = getTimeUntilNextReaction(getTotalPropensity());
            if(currentTime + timeUntilNextReaction > endValue && !endInfinity)
                break;
            updateTime(timeUntilNextReaction);
            allTimes.push_back(currentTime);
        }
        double uniformRV = getUniformRandomVariable();
        int reactionIndex = reaction_tree->searchForNode(uniformRV);
        updateState(stateChangeVector, reactionIndex);
        if (!finalOnly) {
            allStates.push_back(currentState);
        }
        vector<int> dependentReactionIndices = dependency_graph->getDependentReactions(reactionIndex);
        for(int reaction: dependentReactionIndices){
            reaction_tree->updatePropensity(reaction, reactionRates[reaction], currentState, reactantsVector[reaction]);
        }
        currentIteration++;
    }
    // If we pass endValue, remove the last time and state before sending
}
/* backend end */

/* Return the version of Library Link */
DLLEXPORT mint WolframLibrary_getVersion() { return WolframLibraryVersion; }

/* Initialize Library */
DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData) { return LIBRARY_NO_ERROR; }

/* Uninitialize Library */
DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData) { return; }

/* convert a 2-dimentional NumericArray to a 2-dimentional vector */
template <typename Tin, typename Tout>
static vector<vector<Tout> > numericMatrixtoVector(const void *in0, mint const *dims) {
	const Tin *in = static_cast<const Tin *>(in0);
	vector<vector<Tout> > out;
	mint row = dims[0];
	mint col = dims[1];
	for (mint i = 0; i < row; i++) {
		vector<Tout> out_row;
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
        out.push_back((Tout)in[i]);
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

/* convert a matrix to a 1-dimentional NumericArray */
template <typename Tin, typename Tout>
static void matrixtoNumericArray(void *Mout0, vector<vector<Tin> > out) {
	Tout *Mout = static_cast<Tout *>(Mout0);
	int64_t row = out.size();
	int64_t col = out[0].size();
	for (mint i = 0; i < row; i++) {
		for (mint j = 0; j < col; j++) {
			Mout[i*col + j] = (Tout)out[i][j];
		}
	}
}

/* construct reactants and state_change array */
template <typename T1, typename T2>
static void reactantsAndStateChangeArrayConstruction(mint reactionCount, mint moleculeCount, const int64_t *reactIn, const int64_t *prodIn, vector<vector<pair<T1, T2> > >& reactantsArray, vector<vector<pair<T1, T2> > >& stateChangeArray) {
	for (mint i = 0; i < reactionCount; i++) {
		vector<pair<T1, T2> > reactantsArray_row;
		vector<pair<T1, T2> > stateChangeArray_row;
		for (mint j = 0; j < moleculeCount; j++) {
			T1 index = (T1)j;
			T2 in = (T2)reactIn[i*moleculeCount + j];
			T2 out = (T2)prodIn[i*moleculeCount + j];
			if(in > 0) {
				reactantsArray_row.push_back(pair<T1, T2>(index, in));
				stateChangeArray_row.push_back(pair<T1, T2>(index, -in));
			}
			if(out > 0) {
                if (in > 0) {
				    stateChangeArray_row.back().second += out;
                } else {
                    stateChangeArray_row.push_back(pair<T1, T2>(index, out));
                }
			}
		}
        reactantsArray.push_back(reactantsArray_row);
		stateChangeArray.push_back(stateChangeArray_row);
    }
}

// ****** gloabl storage for returns *******
vector<vector<int> > allStates;
vector<double> allTimes;

// ******** end of global storage ********

/* CRN SSA main function */
EXTERN_C DLLEXPORT int directSSAInterface(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument res) {
	// debug setup
	int err = LIBRARY_FUNCTION_ERROR;
	WolframNumericArrayLibrary_Functions naFuns = libData->numericarrayLibraryFunctions;

	// convert initCounts
	MNumericArray MinitCounts = MArgument_getMNumericArray(Args[0]);
	void* MinitCounts_in = naFuns->MNumericArray_getData(MinitCounts);
	mint length = naFuns->MNumericArray_getFlattenedLength(MinitCounts);
    const int64_t *initIn = static_cast<const int64_t *>(MinitCounts_in);
	vector<int> moleculeAmounts = numericArraytoVector<int, int>(initIn, length);

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
    vector<vector<pair<int, int> > > reactantsArray;
	vector<vector<pair<int, int> > > stateChangeArray;
    reactantsAndStateChangeArrayConstruction<int, int>(reactionCount, moleculeCount, reactIn, prodIn, reactantsArray, stateChangeArray);

	// convert rates
	MNumericArray Mrates = MArgument_getMNumericArray(Args[3]);
	void* rateIn = naFuns->MNumericArray_getData(Mrates);
	length = naFuns->MNumericArray_getFlattenedLength(Mrates);
	vector<double> kValues = numericArraytoVector<double, double>(rateIn, length);

	// extract timeEndR
	double timeEndR = MArgument_getReal(Args[4]);

    // extract iterEndI
    double iterEndI = (double)MArgument_getInteger(Args[5]);

    // extract inf (flag)
    mbool inf = MArgument_getBoolean(Args[6]);

    // extract useIter (flag)
    mbool useIter = MArgument_getBoolean(Args[7]);

    // extract statesOnly (flag)
    mbool statesOnly = MArgument_getBoolean(Args[8]);

    // extract finalOnly (flag)
    mbool finalOnly = MArgument_getBoolean(Args[9]);

    // choose endValue
    double endValue = (useIter)? iterEndI : timeEndR;

	// CRN SSA process: pass everything to backend
	directMethodSSA* process = new directMethodSSA(
					moleculeAmounts,
					kValues,
					reactantsArray,
					stateChangeArray,
					endValue,
                    statesOnly,
                    finalOnly,
                    inf,
                    useIter);
	process->start();

    // pass back rerults depends on result
    if (finalOnly) {
		vector<vector<int> > current_state;
		current_state.push_back(process->getCurrentState());
        allStates = current_state;
        if (!statesOnly) {
			vector<double> current_time;
			current_time.push_back(process->getCurrentTime());
            allTimes = current_time;
        }
    } else {
        allStates = process->getAllStates();
        if (!statesOnly) {
            allTimes = process->getAllTimes();
        }
    }

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
	const mint *dims_out = out_size;
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
