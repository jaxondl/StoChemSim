/* Include required headers */
#include <cstdint>
#include <vector>

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
static void vectortoNumericArray(void *Mout0, std::vector<double> out) {
	double *Mout =  static_cast<double *>(Mout0);
	for (int64_t i = 0; i < out.size(); i++) {
		Mout[i] = out[i];
	}
}

/* test function */
static std::vector<double> test(std::vector<std::int64_t> array, std::vector<std::vector<std::int64_t>> matrix, int64_t constant) {
    for(int i = 0; i < matrix.size(); i++) {
        array.insert(array.end(), matrix[i].begin(), matrix[i].end());
    }
    std::vector<double> ans;
    for(int i = 0; i < array.size(); i++) {
        ans.push_back(array[i]/(double)constant);
    }
    return ans;
}

/* CRN SSA main function */
EXTERN_C DLLEXPORT int testing(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument res) {
	// debug setup
	int err = LIBRARY_FUNCTION_ERROR;
	WolframNumericArrayLibrary_Functions naFuns = libData->numericarrayLibraryFunctions;

	// reused local varibles setup
	void *data_in = NULL, *data_out = NULL;
	mint length;
	mint const *dims;

	// matrix conversion test
	MNumericArray Marray = MArgument_getMNumericArray(Args[0]);
	data_in = naFuns->MNumericArray_getData(Marray);
	length = naFuns->MNumericArray_getFlattenedLength(Marray);
	std::vector<std::int64_t> array = numericArraytoVector(data_in, length);

	// matrix conversion test
	MNumericArray Mmatrix = MArgument_getMNumericArray(Args[1]);
	data_in = naFuns->MNumericArray_getData(Mmatrix);
	dims = naFuns->MNumericArray_getDimensions(Mmatrix);
	std::vector<std::vector<std::int64_t>> matrix = numericMatrixtoVector(data_in, dims);

	// convert tEnd
	mint Mnumber = MArgument_getInteger(Args[2]);
	int64_t number = (int64_t)Mnumber;

	// CRN SSA process: pass everything to backend
    std::vector<double> out = test(array, matrix, number);

	// output setup
	MNumericArray Mout;
	int64_t out_size = out.size();
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
	vectortoNumericArray(data_out, out);

	// pass the result back
	MArgument_setMNumericArray(res, Mout);
	return LIBRARY_NO_ERROR;

	cleanup:
	naFuns->MNumericArray_free(Mout);
	return err;
}
