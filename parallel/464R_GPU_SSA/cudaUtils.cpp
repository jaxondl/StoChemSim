#include "cudaUtils.h"

void save_config(std::string fname, std::vector<int> config) {
	std::ofstream out(fname, std::ios::binary);
	uint64_t size = config.size();
	out.write(reinterpret_cast<char*>(&size), sizeof(size));
	out.write(reinterpret_cast<char*>(config.data()), size * sizeof(int));
	out.close();
}

void save_times(std::string fname, std::vector<double> times) {
	std::ofstream out(fname, std::ios::binary);
	uint64_t size = times.size();
	out.write(reinterpret_cast<char*>(&size), sizeof(size));
	out.write(reinterpret_cast<char*>(times.data()), size * sizeof(double));
	out.close();
}

bool is_stable(bool* stability_flags, int n) {
	for (int idx = 0; idx < n; idx++) {
		if (stability_flags[idx] != true) {
			return false;
		}
	}
	return true;
}