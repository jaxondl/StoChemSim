/* Include required headers */
#include <cstdint>
#include <vector>

#include "WolframLibrary.h"
#include "WolframNumericArrayLibrary.h"

/* function prototypes for the backend team 
   TODO: make it a header file, and link it
*/

void preprocess_input (	std::vector<int64_t> initCounts,
						std::vector<std::vector<int64_t>> reactCounts, 
						std::vector<std::vector<int64_t>> prodCounts, 
						std::vector<int64_t> rates, 
						int64_t tEnd);
std::vector<int64_t> run_SSA();

/* =========backend code start=======
*  please follow the fucntion prototype
*  =========backend code end=======/

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
static std::vector<std::vector<int64_t>> numericMatrixtoVector(const void *in0, mint const *dims) {
	const int64_t *in = static_cast<const int64_t *>(in0);
	std::vector<std::vector<int64_t>> out;
	mint row = dims[0];
	mint col = dims[1];
	for (mint i = 0; i < row; i++) {
		std::vector<int64_t> out_row;
		for (mint j = 0; j < col; j++) {
			out_row.push_back(in[i*col + j]);
		}
		out.push_back(out_row);
	}
	return out;
}

/* convert a 1-dimentional NumericArray to a 1-dimentional vector */
static std::vector<int64_t> numericArraytoVector(const void *in0, mint const length) {
	const int64_t *in = static_cast<const int64_t *>(in0);
    std::vector<int64_t> out;
    for (mint i = 0; i < length; i++) {
        out.push_back(in[i]);
    }
	return out;
}

/* convert a 1-dimentional vector to a 1-dimentional NumericArray */
static void vectortoNumericArray(void *Mout0, std::vector<int64_t> out) {
	int64_t *Mout =  static_cast<int64_t *>(Mout0);
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
	void *data_in = NULL, *data_out = NULL;
	mint length;
	mint const *dims;

	// convert initCounts
	MNumericArray MinitCounts = MArgument_getMNumericArray(Args[0]);
	data_in = naFuns->MNumericArray_getData(MinitCounts);
	length = naFuns->MNumericArray_getFlattenedLength(MinitCounts);
	std::vector<int64_t> initCounts = numericArraytoVector(data_in, length);

	// convert reactCounts
	MNumericArray MreactCounts = MArgument_getMNumericArray(Args[1]);
	data_in = naFuns->MNumericArray_getData(MreactCounts);
	dims = naFuns->MNumericArray_getDimensions(MreactCounts);
	std::vector<std::vector<int64_t>> reactCounts = numericMatrixtoVector(data_in, dims);

	// convert prodCounts
	MNumericArray MprodCounts = MArgument_getMNumericArray(Args[2]);
	data_in = naFuns->MNumericArray_getData(MprodCounts);
	dims = naFuns->MNumericArray_getDimensions(MprodCounts);
	std::vector<std::vector<int64_t>> prodCounts = numericMatrixtoVector(data_in, dims);

	// convert rates
	MNumericArray Mrates = MArgument_getMNumericArray(Args[3]);
	data_in = naFuns->MNumericArray_getData(Mrates);
	length = naFuns->MNumericArray_getFlattenedLength(Mrates);
	std::vector<std::int64_t> rates = numericArraytoVector(data_in, length);

	// convert tEnd
	mint MtEnd = MArgument_getInteger(Args[4]);
	int64_t tEnd = (int64_t)MtEnd;

	// CRN SSA process: pass everything to backend
	preprocess_input(initCounts, reactCounts, prodCounts, rates, tEnd);
	std::vector<int64_t> out = run_SSA();

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
