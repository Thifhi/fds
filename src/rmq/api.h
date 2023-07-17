#include <vector>
#include <cstdint>
#include <chrono>

#include "rmq.h"

using namespace std;

// Naive implementation
pair<vector<uint64_t>, uint64_t> naive(vector<uint64_t>& data, vector<pair<uint64_t, uint64_t>>& queries) {
    auto index = in_block_rmq(data);

    vector<uint64_t> results(queries.size());
    for (uint64_t i = 0; i < queries.size(); ++i) {
        auto& query = queries[i];
        auto res = index[get_index_from_s_e(query.first, query.second, data.size())];
        results[i] = res;
    }
    return make_pair(results, 0);

}

// N-Logn implementation
pair<vector<uint64_t>, uint64_t> nlogn_rmq(vector<uint64_t>& data, vector<pair<uint64_t, uint64_t>>& queries) {
    auto index_structure = generate_power_of_two_index_struct(data);
    
    vector<uint64_t> results(queries.size());
    for (uint64_t i = 0; i < queries.size(); ++i) {
        auto& query = queries[i];
        auto res = query_blockwise(query.first, query.second, data, index_structure);
        results[i] = res;
    }
    return make_pair(results, 0);
}

// Linear implementation
pair<vector<uint64_t>, uint64_t> linear_rmq(vector<uint64_t>& data, vector<pair<uint64_t, uint64_t>>& queries) {
    // Keep track of the used space to be reported
    uint64_t bitsize = 0;

    uint64_t block_size = ceil(log2(data.size()) / 4.0);

    // Append UINT64_MAXs to avoid partial block in the end
    auto num_blocks = ceil(static_cast<double>(data.size()) / block_size);
    data.reserve(block_size * num_blocks);
    for (uint64_t i = data.size(); i < block_size * num_blocks; ++i) {
        data.emplace_back(UINT64_MAX);
    }

    // Build the cartesian tree per block
    vector<vector<bool>> block_cartesian_trees(num_blocks);
    for (uint64_t i = 0; i < num_blocks; ++i) {
        vector<uint64_t> slice = vector<uint64_t>(data.begin() + i * block_size, data.begin() + (i + 1) * block_size);
        auto block_tree = build_cartesian_tree(slice);
        block_cartesian_trees[i] = block_tree;
    }
    for (auto& item : block_cartesian_trees) {
        bitsize += item.size();
    }

    // Store the minimum of each block for block-wise queries
    vector<uint64_t> block_mins(num_blocks);
    // Store the actual indices of the minimums so that we can map back
    vector<uint64_t> block_min_indices(num_blocks);
    for (uint64_t i = 0; i < num_blocks; ++i) {
        uint64_t min = UINT64_MAX;
        uint64_t min_pos;
        for (uint64_t j = 0; j < block_size; ++j) {
            if (data[i * block_size + j] < min) {
                min = data[i * block_size + j];
                min_pos = i * block_size + j;
            }
        }
        block_min_indices[i] = min_pos;
        block_mins[i] = min;
    }

    bitsize += block_min_indices.size() * 64;
    bitsize += block_min_indices.size() * 64;

    auto lookup_table = generate_lookup_table(block_size);
    for (auto& items : lookup_table) {
        bitsize += items.second.size() * 64;
    }

    auto index_struct = generate_power_of_two_index_struct(block_mins);
    bitsize += index_struct.first.size() * 64;
    
    // Handle queries one-by-one
    vector<uint64_t> results(queries.size());
    for (uint64_t i = 0; i < queries.size(); ++i) {
        auto s = queries[i].first;
        auto e = queries[i].second;
        uint64_t s_block_idx = s / block_size;
        uint64_t e_block_idx = e / block_size;
        uint64_t sl = s % block_size;
        uint64_t er = e % block_size;
        uint64_t res;
        if (s_block_idx == e_block_idx) {
            // The query s-e is both in the same block
            res = query_inside_block(sl, er, block_cartesian_trees, s_block_idx, block_size, lookup_table);
            res += s_block_idx * block_size;
        } else {
            vector<uint64_t> possible_indices;
            // Query for the block that s belongs in
            uint64_t left = query_inside_block(sl, block_size - 1, block_cartesian_trees, s_block_idx, block_size, lookup_table);
            left += s_block_idx * block_size;
            possible_indices.emplace_back(left);
            // Query for the block that e belongs in
            uint64_t right = query_inside_block(0, er, block_cartesian_trees, e_block_idx, block_size, lookup_table);
            right += e_block_idx * block_size;
            possible_indices.emplace_back(right);
            // Query for the block in-between (if there are any)
            if (e_block_idx - s_block_idx > 1) {
                uint64_t block_idx = query_blockwise(s_block_idx + 1, e_block_idx - 1, block_mins, index_struct);
                possible_indices.emplace_back(block_min_indices[block_idx]);
            }
            uint64_t min = UINT64_MAX;
            uint64_t min_idx;
            for (uint64_t i = 0; i < possible_indices.size(); ++i) {
                if (data[possible_indices[i]] < min) {
                    min = data[possible_indices[i]];
                    min_idx = possible_indices[i];
                }
            }
            res = min_idx;
        }
        results[i] = res;
    }
    return make_pair(results, bitsize);
}

auto solve_rmq(vector<uint64_t>& data, vector<pair<uint64_t, uint64_t>>& queries, string algorithm_type = "linear") {
    auto start = std::chrono::high_resolution_clock::now();
    pair<vector<uint64_t>, uint64_t> res;

    if (algorithm_type == "naive") {
        res = naive(data, queries);
    } else if (algorithm_type == "nlogn") {
        res = nlogn_rmq(data, queries);
    } else if (algorithm_type == "linear") {
        res = linear_rmq(data, queries);
    } else {
        cout << "Unknown algorithm!" << endl;
        exit(1);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    return make_tuple(res.first, duration, res.second);
}