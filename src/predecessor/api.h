#include <vector>
#include <cstdint>
#include <fstream>
#include "predecessor.h"

using namespace std;

vector<uint64_t> solve_pd(vector<uint64_t>& data, vector<uint64_t>& queries) {
    auto index_structure = EliasFanoCoded(data);

    vector<uint64_t> results(queries.size());
    for (auto query : queries) {
        auto idx = index_structure.predecessor_idx(query);
        results.emplace_back(data[idx]);
    }
    return results;
}