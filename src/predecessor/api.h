#include <vector>
#include <cstdint>
#include <chrono>

#include "predecessor.h"

using namespace std;

auto solve_pd(vector<uint64_t>& data, vector<uint64_t>& queries) {
    auto start = chrono::high_resolution_clock::now();

    auto index_structure = EliasFanoCoded(data);
    vector<uint64_t> results(queries.size());
    for (uint64_t i = 0; i < queries.size(); ++i) {
        results[i] = index_structure.predecessor(queries[i]);
    }

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();

    uint64_t bits = index_structure.get_size();

    return make_tuple(results, duration, bits);
}