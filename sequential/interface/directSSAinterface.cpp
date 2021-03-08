/* Include required headers */
#include <cstdint>
#include <vector>
#include "../directSSA/directMethodSSA.h"
#include "WolframLibrary.h"
#include "WolframNumericArrayLibrary.h"

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
static vector<vector<Tout>> numericMatrixtoVector(const void *in0, mint const *dims) {
	const Tin *in = static_cast<const Tin *>(in0);
	vector<vector<Tout>> out;
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

template <typename Tin, typename Tout>
static void matrixtoNumericArray(void *Mout0, vector<vector<Tin>> out) {
	Tout *Mout = static_cast<Tout *>(Mout0);
	int64_t row = out.size();
	int64_t col = out[0].size();
	for (mint i = 0; i < row; i++) {
		for (mint j = 0; j < col; j++) {
			Mout[i*col + j] = (Tout)out[i][j];
		}
	}
}

template <typename T1, typename T2>
static void reactantsAndStateChangeArrayConstruction(mint reactionCount, mint moleculeCount, const int64_t *reactIn, const int64_t *prodIn, vector<vector<pair<T1, T2>>>& reactantsArray, vector<vector<pair<T1, T2>>>& stateChangeArray) {
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
vector<vector<int>> allStates;
vector<double> allTimes;

// ******** end of global storage ********

/* CRN SSA main function */
EXTERN_C DLLEXPORT int CRN_SSA(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument res) {
	// debug setup
	int err = LIBRARY_FUNCTION_ERROR;
	WolframNumericArrayLibrary_Functions naFuns = libData->numericarrayLibraryFunctions;

	// convert initCounts
	MNumericArray MinitCounts = MArgument_getMNumericArray(Args[0]);
	void* MinitCounts_in = naFuns->MNumericArray_getData(MinitCounts);
	mint length = naFuns->MNumericArray_getFlattenedLength(MinitCounts);
    const int64_t *initIn = static_cast<const int64_t *>(MinitCounts_in);   // TODO: convert to template in later release
	vector<int> moleculeAmounts = numericArraytoVector<int, int>(initIn, length);

	// convert reactantsArray & stateChangeArray
	MNumericArray MreactCounts = MArgument_getMNumericArray(Args[1]);
	MNumericArray MprodCounts = MArgument_getMNumericArray(Args[2]);
	mint const * dims = naFuns->MNumericArray_getDimensions(MreactCounts);
	mint reactionCount = dims[0];
	mint moleculeCount = dims[1];
	void* MreactCounts_in = naFuns->MNumericArray_getData(MreactCounts);
	void* MprodCounts_in = naFuns->MNumericArray_getData(MprodCounts);
	const int64_t *reactIn = static_cast<const int64_t *>(MreactCounts_in); // TODO: convert to template in later release
	const int64_t *prodIn = static_cast<const int64_t *>(MprodCounts_in);   // TODO: convert to template in later release
    vector<vector<pair<int, int>>> reactantsArray;
	vector<vector<pair<int, int>>> stateChangeArray;
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

    if (finalOnly) {
		vector<vector<int>> current_state;
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
