#include <vector>
#include <cstdint>

#include "predecessor.h"

using namespace std;

vector<uint64_t> solve_pd(vector<uint64_t>& data, vector<uint64_t>& queries) {
    auto index_structure = EliasFanoCoded(data);

    vector<uint64_t> results(queries.size());
    for (uint64_t i = 0; i < queries.size(); ++i) {
        auto idx = index_structure.predecessor_idx(queries[i]);
        if (idx == UINT64_MAX) {
            results[i] = UINT64_MAX;
        } else {
            results[i] = (data[idx]);
        }
    }
    return results;
}