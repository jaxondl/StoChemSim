/* Include required headers */
#include <cstdint>
#include <vector>

#include "WolframLibrary.h"
#include "WolframNumericArrayLibrary.h"

/* function prototypes for the backend team 
   TODO: make it a header file, and link it
   COMMENT: I put uint32_t instead of uint8_t in the second parameter of pair for now
*/

void preprocess_input (	std::vector<uint32_t> moleculeAmounts,
						std::vector<std::vector<std::pair<uint32_t, uint32_t>>> stateChangeArray,
						std::vector<std::vector<std::pair<uint32_t, uint32_t>>> reactantsArray,
						std::vector<double> kValues,
						double tEnd);
std::vector<double> run_SSA();

// =========backend code start=======

std::vector<double> run_SSA() {
	std::vector<double> ans;
	ans.push_back(1.2);
	ans.push_back(2.4);
	ans.push_back(1);
	return ans;
}

// =========backend code end=======

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
static std::vector<std::vector<int64_t>> numericMatrixtoVector(const void *in0, mint const *dims) {
	const Tin *in = static_cast<const Tin *>(in0);
	std::vector<std::vector<Tout>> out;
	mint row = dims[0];
	mint col = dims[1];
	for (mint i = 0; i < row; i++) {
		std::vector<T> out_row;
		for (mint j = 0; j < col; j++) {
			out_row.push_back(in[i*col + j]);
		}
		out.push_back(out_row);
	}
	return out;
}

/* convert a 1-dimentional NumericArray to a 1-dimentional vector */
template <typename Tin, typename Tout>
static std::vector<Tout> numericArraytoVector(const void *in0, mint const length) {
	const Tin *in = static_cast<const Tin *>(in0);
    std::vector<Tout> out;
    for (mint i = 0; i < length; i++) {
        out.push_back(in[i]);
    }
	return out;
}

/* convert a 1-dimentional vector to a 1-dimentional NumericArray */
template <typename T>
static void vectortoNumericArray(void *Mout0, std::vector<T> out) {
	T *Mout =  static_cast<T *>(Mout0);
	for (int64_t i = 0; i < out.size(); i++) {
		Mout[i] = out[i];
	}
}

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
	std::vector<uint32_t> moleculeAmounts = numericArraytoVector<int64_t, uint32_t>(data_in, length);


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

    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> reactantsArray;
	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> stateChangeArray;
    for (mint i = 0; i < reactionCount; i++) {
		std::vector<std::pair<uint32_t, uint32_t>> reactantsArray_row;
		std::vector<std::pair<uint32_t, uint32_t>> stateChangeArray_row;
		for (mint j = 0; j < moleculeCount; j++) {
			uint32_t index = (uint32_t)j;
			uint32_t in = (uint32_t)reactIn[i*moleculeCount + j];
			uint32_t out = (uint32_t)prodIn[i*moleculeCount + j];
			if(in > 0) {
				reactantsArray_row.push_back(std::pair<uint32_t, uint32_t>(index, in));
				stateChangeArray_row.push_back(std::pair<uint32_t, uint32_t>(index, -in));
			}
			if(out > 0) {
				stateChangeArray_row.push_back(std::pair<uint32_t, uint32_t>(index, out));
			}
		}
        reactantsArray.push_back(reactantsArray_row);
		stateChangeArray.push_back(stateChangeArray_row);
    }

	// convert rates
	MNumericArray Mrates = MArgument_getMNumericArray(Args[3]);
	data_in = naFuns->MNumericArray_getData(Mrates);
	length = naFuns->MNumericArray_getFlattenedLength(Mrates);
	std::vector<double> kValues = numericArraytoVector<double, double>(data_in, length);

	// convert tEnd
	mreal tEnd = MArgument_getReal(Args[4]);

	// CRN SSA process: pass everything to backend
	preprocess_input(moleculeAmounts, stateChangeArray, reactantsArray, kValues, tEnd);
	std::vector<double> out = run_SSA();

	// output setup
	MNumericArray Mout;
	int64_t out_size = out.size();
	const mint *dims_out = &out_size;
	err = naFuns->MNumericArray_new(MNumericArray_Type_Bit64, 1, dims_out, &Mout);
	if (err != 0) {
		goto cleanup;
	}
	data_out = naFuns->MNumericArray_getData(Mout);
	if (data_out == NULL) {
		goto cleanup;
	}
	
	// convert output to a NumericArray
	vectortoNumericArray(data_out, out);

	// pass the result back
	MArgument_setMNumericArray(res, Mout);
	return LIBRARY_NO_ERROR;

	cleanup:
	naFuns->MNumericArray_free(Mout);
	return err;
}
