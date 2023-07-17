#include <memory>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <algorithm>
#include <cmath>

using namespace std;

// Builds the per-block cartesian tree using Cartesian Tree Signature encoding from:
// https://link.springer.com/article/10.1007/s00453-012-9683-x
vector<bool> build_cartesian_tree(vector<uint64_t>& block) {
    vector<bool> succ_cartesian_tree(2 * block.size());
    uint64_t tree_pointer = 0;
    vector<uint64_t> stack(block.size());
    uint64_t stack_pointer = UINT64_MAX;
    for (int i = 0; i < block.size(); ++i) {
        while (stack_pointer != UINT64_MAX && stack[stack_pointer] > block[i]) {
            succ_cartesian_tree[tree_pointer++] = 1;
            --stack_pointer;
        }
        succ_cartesian_tree[tree_pointer++] = 0;
        stack[++stack_pointer] = block[i];
    }
    succ_cartesian_tree.resize(tree_pointer);
    return succ_cartesian_tree;
}

// Map the start - end indices into a single value so that we can use flat vectors instead of nested
uint64_t get_index_from_s_e(uint64_t s, uint64_t e, uint64_t block_size) {
    uint64_t n = block_size - (s - 1);
    uint64_t m = block_size;
    uint64_t idx = (e - s) + ((n + m) / 2.0) * (m - n + 1);
    return idx;
}

// Naive RMQ that is used for per block queries, or the whole data for a naive implementation
// Returned data structure contains the RMQ for each (start, end) combination in a flat-vector
vector<uint64_t> in_block_rmq(vector<uint64_t>& block) {
    // Block size is tiny with log(n) / 4, so naive implementation should be fast
    vector<uint64_t> res(block.size() * (block.size() + 1) / 2);
    for (uint64_t i = 0; i < block.size(); ++i) {
        uint64_t min = UINT64_MAX;
        uint64_t min_pos;
        for (uint64_t j = i; j < block.size(); ++j) {
            if (block[j] < min) {
                min = block[j];
                min_pos = j;
            }
            auto idx = get_index_from_s_e(i, j, block.size());
            res[idx] = min_pos;
        }
    }
    return res;
}

// Generation of the lookup data structure for per block RMQ with (Tree, start, end) triplet
unordered_map<vector<bool>, vector<uint64_t>> generate_lookup_table(uint64_t block_size) {
    unordered_map<vector<bool>, vector<uint64_t>> table;
    vector<uint64_t> permutate(block_size); 

    for (uint64_t i = 0; i < permutate.size(); ++i) {
        permutate[i] = i;
    }

    // For each possible block permutation
    do {
        auto tree = build_cartesian_tree(permutate);
        // Tree is not unique, so skip if we already added the same tree before
        if (table.find(tree) != table.end()) {
            continue;
        }
        auto permutation_rmq = in_block_rmq(permutate);
        table[tree] = permutation_rmq;
    } while (next_permutation(permutate.begin(), permutate.end()));
    return table;
}

uint64_t query_inside_block(uint64_t s, uint64_t e, vector<vector<bool>>& block_trees, uint64_t block_idx, uint64_t block_size, unordered_map<vector<bool>, vector<uint64_t>>& lookup_table) {
    if (s == e) {
        return s;
    }
    auto& our_block_tree = block_trees[block_idx];
    return lookup_table[our_block_tree][get_index_from_s_e(s, e, block_size)];
}

// Build the nlogn space index data structure
pair<vector<uint64_t>, uint64_t> generate_power_of_two_index_struct(vector<uint64_t>& data) {
    // Max look-ahead stored is 2^k
    uint64_t k = ceil(log2(data.size()));
    // Store the lookup table in a flat vector and access i-th element with i / k + i % k
    vector<uint64_t> power_of_two_index_struct(k * data.size());

    // Fill the first entry
    for (uint64_t i = 0; i < data.size() - 1; ++i) {
        power_of_two_index_struct[i * k] = (data[i] < data[i + 1]) ? i : i + 1;
    }
    power_of_two_index_struct[(data.size() - 1) * k] = data.size() - 1;

    // Fill the rest via dynamic programming
    for (uint64_t j = 1; j < k; ++j) {
        for (uint64_t i = 0; i < data.size(); ++i) {
            if (i + pow(2, j) >= data.size()) {
                power_of_two_index_struct[i * k + j] = power_of_two_index_struct[i * k + j - 1];
            } else {
                uint64_t right_half_min_idx = power_of_two_index_struct[(i + pow(2, j)) * k + j - 1];
                uint64_t left_half_min_idx = power_of_two_index_struct[i * k + j - 1];
                if (data[left_half_min_idx] < data[right_half_min_idx]) {
                    power_of_two_index_struct[i * k + j] = left_half_min_idx;
                } else {
                    power_of_two_index_struct[i * k + j] = right_half_min_idx;
                }
            }
        }
    }
    return make_pair(power_of_two_index_struct, k);
}

// Return the arg-min of rmq(s, 2^l) and rmq(e - 2^l + 1, e)
uint64_t query_blockwise(uint64_t s, uint64_t e, vector<uint64_t>& data, pair<vector<uint64_t>, uint64_t>& index_struct) {
    if (e - s <= 1) {
        return (data[s] < data[e]) ? s : e;
    }
    auto& index = index_struct.first;
    auto k = index_struct.second;
    uint64_t l = floor(log2(e - s));
    auto first = index[s * k + l - 1];
    auto second = index[(e - pow(2, l) + 1) * k + l - 1];
    return (data[first] < data[second]) ? first : second;
}
