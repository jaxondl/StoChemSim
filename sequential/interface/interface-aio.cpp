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
    vector<vector<int> > dependencyGraphStructure; // the graph structure itself
public:
    dependencyGraph(vector<vector<pair<int, int> > > stateChangeVector, vector<vector<pair<int,int> > > reactantsVector, vector<int> moleculeAmounts); // constructor

    bool intersects(set<int> set1, set<int> set2); // user made function to check for an intersection of 1 or more elements between two sets

    vector<int> getDependentReactions(int reactionIndex); // obtain the dependent reactions for a given reaction via the dependencyGraphStructure
};

dependencyGraph::dependencyGraph(vector<vector<pair<int, int> > > stateChangeVector, vector<vector<pair<int,int> > > reactantsVector, vector<int> moleculeAmounts){
    // DEFINITION 1: Reactants(p) and products(p) are reactants and prods of reaction p. e.g. Reactants(1) = {a,b}
    // DEFINITION 2: DependsOn(a-mu), where a-mu is the propensity of chosen reaction, is the set of substances that affect its value. i.e. Reactants(mu)
    // DEFINITION 3: Affects(mu) is set of substances that change in number when reaction mu executes. i.e. Reactants(mu) UNION Products(mu)
    /** DEFINITION 4 (DEPENDENCY GRAPH): Directed graph where vertex set=R (all reactions),
     * directed edge FROM Vi to Vj IFF Affects(Vi) INTERSECTION DependsOn(a-Vj) is NOT an empty set.
     * i.e. at least one of the reactants and products of Vi is shared with the reactants of Vj.
    */

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

// Propesnity calculation for each reaction = k * (amount(reactantMolecule1)) * (amount(reactantMolecule1)-1)...(amount(reactantMolecule1)-reactantCoefficient + 1) * ... continue for every reactant
double reactionTree::calculatePropensity(double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int> > reactants){
    double propensity = reactionRate;
    for (pair<int, int> reactant: reactants) { 
        for (int i=0; i<reactant.second; i++) {
            propensity *= (moleculeAmounts[reactant.first] - i);
        }
    }
    return propensity;
}

// Create the reaction tree
reactionTree::reactionTree(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector) {
    int numReactions = reactantsVector.size();
    reactionTreeArray = new reactionNode[numReactions];
    reactionTreeArray[0].parent = -1;
    // Calculate propensity for each reaction
    for (int i=0; i<numReactions; i++) {
        reactionTreeArray[i].propensity = calculatePropensity(reactionRates[i], moleculeAmounts, reactantsVector[i]);
        reactionTreeArray[i].leftSum = 0;
        reactionTreeArray[i].rightSum = 0;
    }
    // Assign left and right child indices for each tree node beginning at the roott
    for (int i=0; i<numReactions; i++) {
        if (i * 2 + 1 < numReactions) { 
            reactionTreeArray[i].leftChild = i * 2 + 1;
            reactionTreeArray[i * 2 + 1].parent = i;
        }
        else {
            reactionTreeArray[i].leftChild = -1;
        }
        if (i * 2 + 2 < numReactions) {
            reactionTreeArray[i].rightChild = i * 2 + 2;
            reactionTreeArray[i * 2 + 2].parent = i;
        }
        else {
            reactionTreeArray[i].rightChild = -1;
        }
    }

    // Calculating left and right sums for all tree nodes starting from the leaves (end of array to the beginning)
    for (int i=numReactions-1; i>0; i--) {
        double subTotalPropensity = reactionTreeArray[i].propensity + reactionTreeArray[i].rightSum + reactionTreeArray[i].leftSum; // current node's subtree propensity sum
        // Add subTotalPropensity to parent's leftSum or rightSum depending on if the child is the left or right child of its parent
        if (i == reactionTreeArray[reactionTreeArray[i].parent].leftChild) {
            reactionTreeArray[reactionTreeArray[i].parent].leftSum += subTotalPropensity;
        }
        else {
            reactionTreeArray[reactionTreeArray[i].parent].rightSum += subTotalPropensity;
        }
    }
}

// Search for node in reaction tree based on sampled uniform RV
int reactionTree::searchForNode(double RV) {
    int currentIndex = 0; // index of checkNode in the reactionTree array, start with root node
    reactionNode checkNode = reactionTreeArray[0]; // node object of the reaction to be checked
    double leftSumTotal = reactionTreeArray[0].leftSum;
    double totalPropensity = reactionTreeArray[0].leftSum + reactionTreeArray[0].rightSum + reactionTreeArray[0].propensity;
    while(RV < (leftSumTotal/totalPropensity) || RV > ((leftSumTotal+checkNode.propensity)/totalPropensity)) { // if the sampled RV is not in the current node's propensity range, continue searching
        if (RV < (leftSumTotal/totalPropensity)) { // if the sampled RV is in the left subtree of the current node, update the current node to the left child
            currentIndex = currentIndex*2 + 1;
            leftSumTotal -= checkNode.leftSum;
            if (checkNode.leftChild == -1){
                return -1;
                break;
            }
            checkNode = reactionTreeArray[checkNode.leftChild];
            leftSumTotal += checkNode.leftSum;
        } else { // if the sampled RV is in the right subtree of the current node, update the current node to the right child
            currentIndex = currentIndex*2 + 2;
            leftSumTotal += checkNode.propensity;
            if (checkNode.rightChild == -1){
                return -1;
                break;
            }
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
    double propensityChange = newPropensity - reactionTreeArray[currentIndex].propensity; // determine change in propensity
    reactionTreeArray[currentIndex].propensity = newPropensity;
    while (reactionTreeArray[currentIndex].parent != -1) { // update parents' left or right sums depending on if the child is the left or right child of its parent until you hit the root node
        if (currentIndex == reactionTreeArray[reactionTreeArray[currentIndex].parent].leftChild) {
            reactionTreeArray[reactionTreeArray[currentIndex].parent].leftSum += propensityChange;
        }
        else {
            reactionTreeArray[reactionTreeArray[currentIndex].parent].rightSum += propensityChange;
        }
        currentIndex = reactionTreeArray[currentIndex].parent; // move up the tree to next parent
    }
}

class directMethodSSA {
private:
    reactionTree* reaction_tree; // reaction tree necessary for this simulation algorithm
    dependencyGraph* dependency_graph; // dependency graph necessary for this simulation algorithm

    vector<double> reactionRates; // reaction rates (k) for each reaction
    vector<vector<pair<int, int> > > reactantsVector; // reactants for each reaction
    vector<vector<pair<int, int> > > stateChangeVector; // state change vectors for each reaction
    vector<int> currentState; // the current state / amounts of species in CRN
    vector<vector<int> > allStates; // vector containing all states calculated by the simulation algorithm
    vector<double> allTimes; // vector containing the time it has taken to complete each iteration and all preceding iterations

    double currentTime; // current time tracked for the simulation
    int currentIteration; // current iteration tracked fo the simulation
    double endValue; // end value will either be a ending time or and ending iteration, depending on whether the user sets the endByIteration flag to true. The simulation will stop before it crosses the endValue

    bool endInfinity; // if set to true, the simulation will continue until no more reactions can take place
    bool statesOnly; // if set to true, the simulation will forego any calculations of times (as well as any updating the time) and only prin out the states 
    bool finalOnly; // if set to true, the simulation will only output/print the final iteration of the simulation
    bool endByIteration; // if set to true, the endValue will be used to determine the ending iteration instead of the ending time

    double getUniformRandomVariable(); // obtain a uniform RV dist sample value
    double getTimeUntilNextReaction(double propensity); // obtain the time until the next reaction (in order to update the time)
    void updateTime(double timeUntilNextReaction); // update the existing time
    void updateState(vector<vector<pair<int, int> > > stateChangeVector, int reactionIndex); // update the state of the species
    double getTotalPropensity(); // obtain the total propensity across all reactions
    mt19937 gen_uni;
    mt19937 gen_exp;

public:
    directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool statesOnly, bool finalOnly, bool endInfinity, bool endByIteration); // constructor
    vector<vector<int> > getAllStates(); // obtain vector of all states
    vector<double> getAllTimes(); // obtain vector of all times
    vector<int> getCurrentState(); // obtain vector of the end/most recent state
    double getCurrentTime(); // obtain end/most recent time
    int getCurrentIteration(); // obtain end/most recent iteration
    void start(); // begin CRN SSA simulation
    
};

directMethodSSA::directMethodSSA(vector<int> moleculeAmounts, vector<double> reactionRates, vector<vector<pair<int, int> > > reactantsVector, vector<vector<pair<int, int> > > stateChangeVector, double endValue, bool statesOnly, bool finalOnly, bool endInfinity, bool endByIteration){
    this->reaction_tree = new class reactionTree(moleculeAmounts, reactionRates, reactantsVector); // create reaction tree
    this->dependency_graph = new class dependencyGraph(stateChangeVector, reactantsVector, moleculeAmounts); // create dependency graph
    this->allStates.push_back(moleculeAmounts); // stored initial state in the allStates vector
    this->allTimes.push_back(0); // store the initial time (0) in the allTimes vector
    this->reactionRates = reactionRates; // k reaction constants
    this->reactantsVector = reactantsVector;
    this->stateChangeVector = stateChangeVector;

    this->currentState = moleculeAmounts;
    this->currentTime = 0; // want to set the current time to 0 before beginning
    this->currentIteration = 0; // likewise set current iteration to 0 before beginning

    this->endValue = endValue; // all related booleans/flags are defined by the user
    this->statesOnly = statesOnly;
    this->finalOnly = finalOnly;
    this->endInfinity = endInfinity;
    this->endByIteration = endByIteration;

    random_device rd_uni; // random seed
    mt19937 gen_uni(rd_uni()); // mersenne twister engine
    this->gen_uni = gen_uni;
    random_device rd_exp; // random seed
    mt19937 gen_exp(rd_exp()); 
    this->gen_exp = gen_exp; // mersenne twister engine
    uniform_real_distribution<double> dis(0.0, 1);
    
    if (this->statesOnly && !this->endByIteration) { // if the user requests to only calculate the states but gives a positive finite value for time (i.e. did not set endByIteration to true) then the end time is automatically set to infinity (as finite time will be unapplicable)
        this->endInfinity = true;
    }
    //asdasdasdasdasdas
}

double directMethodSSA::getUniformRandomVariable() {
    //double totalPropensity = getTotalPropensity();
    uniform_real_distribution<double> dis(0.0, 1);
    double RV = dis(gen_uni); // obtain sampled value
    return RV;
}

void directMethodSSA::updateTime(double timeUntilNextReaction){
    currentTime += timeUntilNextReaction; // update time by adding the passed in time until next reaction to the existing time
}

void directMethodSSA::updateState(vector<vector<pair<int, int> > > stateChangeVector, int reactionIndex) {
    vector<pair<int, int> > chosenReactionChange = stateChangeVector[reactionIndex]; // given the chosen reaction, obtain all the changes to all the species
    for(pair<int, int> p: chosenReactionChange){
        currentState[p.first] += p.second; // update the state change by individually updating the amounts of the affected species
    }
}

double directMethodSSA::getTimeUntilNextReaction(double propensity) {
    exponential_distribution<> dis(propensity); // utilize exponential approximation of propensity
    double RV = dis(gen_exp); // obtain sampled value
    return RV;
}

double directMethodSSA::getTotalPropensity(){
    reactionTree::reactionNode root = reaction_tree->reactionTreeArray[0]; // propensity sum = root node propensity + rightSum + leftSum
    double totalPropensity = root.propensity + root.leftSum + root.rightSum;
    return totalPropensity;
}

vector<vector<int> > directMethodSSA::getAllStates(){return allStates;}

vector<double> directMethodSSA::getAllTimes(){return allTimes;}

vector<int> directMethodSSA::getCurrentState(){return currentState;}

double directMethodSSA::getCurrentTime(){return currentTime;}

int directMethodSSA::getCurrentIteration(){return currentIteration;}

void directMethodSSA::start(){
    // continue the simulation while the total propensity > 0 AND (the endInfinity flag is true OR the current time/iteration hasn't exceeded the inputted limit)
    while (getTotalPropensity() > 0 && ((!endByIteration && (currentTime < endValue || endInfinity)) || (endByIteration && (currentIteration < endValue || endInfinity)))){
        if (!statesOnly) { // only calculate if the user wants to also calculate the times (default)
            double timeUntilNextReaction = getTimeUntilNextReaction(getTotalPropensity()); // obtain the time until the next reaction
            if(currentTime + timeUntilNextReaction > endValue && !endInfinity) // if updating the time violates the finite end time value, terminate the simulation
                break;
            updateTime(timeUntilNextReaction); // otherwise, update the time
            allTimes.push_back(currentTime); // recorded the updated time
        }
        double RV = getUniformRandomVariable(); // obtain sampled value
        int reactionIndex = reaction_tree->searchForNode(RV); // search for reaction
        if(reactionIndex == -1)
            break;
        updateState(stateChangeVector, reactionIndex); // update state/configuration
        if (!finalOnly) { // only save the state vectors of the iteration if the finalOnly flag is false
            allStates.push_back(currentState);
        }
        vector<int> dependentReactionIndices = dependency_graph->getDependentReactions(reactionIndex); // get affected reactions from dependency graph
        for(int reaction: dependentReactionIndices){
            reaction_tree->updatePropensity(reaction, reactionRates[reaction], currentState, reactantsVector[reaction]); // update propensity for every affected reaction
        }
        currentIteration++; // update iteration
    }
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
