#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include "rmq.h"

int main_nlogn(int argc, char* argv[]) {
    // std::string inputFileName = "./rmq_example_1.txt";
    std::string inputFileName = argv[1];

    // Open the input file
    std::ifstream input_file(inputFileName);
    if (!input_file) {
        std::cerr << "Failed to open input file: " << inputFileName << std::endl;
        return 1;
    }

    std::string outputFileName = argv[2];
    // Open the output file
    std::ofstream output_file(outputFileName);
    if (!output_file) {
        std::cerr << "Failed to open output file: " << outputFileName << std::endl;
        input_file.close();
        return 1;
    }

    int n;
    input_file >> n;

    vector<uint64_t> data(n);
    for (int i = 0; i < data.size(); ++i) {
        uint64_t num;
        input_file >> num;
        data[i] = num;
    }

    auto index_struct = generate_power_of_two_index_struct(data);

    char comma;
    uint64_t s;
    uint64_t e;
    while (input_file >> s) {
        input_file >> comma;
        input_file >> e;

        auto res = query_blockwise(s, e, data, index_struct);

        output_file << res << std::endl;
    }

    // Close the files
    input_file.close();
    output_file.close();

    return 0;
}

int main(int argc, char* argv[]) {
    // std::string inputFileName = "./rmq_example_1.txt";
    std::string inputFileName = argv[1];

    // Open the input file
    std::ifstream input_file(inputFileName);
    if (!input_file) {
        std::cerr << "Failed to open input file: " << inputFileName << std::endl;
        return 1;
    }

    std::string outputFileName = argv[2];
    // Open the output file
    std::ofstream output_file(outputFileName);
    if (!output_file) {
        std::cerr << "Failed to open output file: " << outputFileName << std::endl;
        input_file.close();
        return 1;
    }

    int n;
    input_file >> n;

    vector<uint64_t> data(n);
    for (int i = 0; i < data.size(); ++i) {
        uint64_t num;
        input_file >> num;
        data[i] = num;
    }

    uint64_t block_size = ceil(log2(data.size()) / 4.0);

    // Avoid partial blocks
    auto num_blocks = ceil(static_cast<double>(data.size()) / block_size);
    data.reserve(block_size * num_blocks);
    for (uint64_t i = data.size(); i < block_size * num_blocks; ++i) {
        data.emplace_back(UINT64_MAX);
    }

    vector<vector<bool>> block_cartesian_trees(num_blocks);
    for (uint64_t i = 0; i < num_blocks; ++i) {
        vector<uint64_t> slice = vector<uint64_t>(data.begin() + i * block_size, data.begin() + (i + 1) * block_size);
        auto block_tree = build_cartesian_tree(slice);
        block_cartesian_trees[i] = block_tree;
    }

    vector<uint64_t> block_min_indices(num_blocks);
    vector<uint64_t> block_mins(num_blocks);
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

    auto lookup_table = generate_lookup_table(block_size);
    auto index_struct = generate_power_of_two_index_struct(block_mins);

    char comma;
    uint64_t s;
    uint64_t e;
    while (input_file >> s) {
        input_file >> comma;
        input_file >> e;
        uint64_t s_block_idx = s / block_size;
        uint64_t e_block_idx = e / block_size;
        uint64_t sl = s % block_size;
        uint64_t er = e % block_size;
        uint64_t res;
        if (s_block_idx == e_block_idx) {
            res = query_inside_block(sl, er, block_cartesian_trees, s_block_idx, block_size, lookup_table);
            res += s_block_idx * block_size;
        } else {
            vector<uint64_t> possible_indices;
            uint64_t left = query_inside_block(sl, block_size - 1, block_cartesian_trees, s_block_idx, block_size, lookup_table);
            left += s_block_idx * block_size;
            possible_indices.emplace_back(left);
            uint64_t right = query_inside_block(0, er, block_cartesian_trees, e_block_idx, block_size, lookup_table);
            right += e_block_idx * block_size;
            possible_indices.emplace_back(right);
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
        output_file << res << std::endl;
    }

    // Close the files
    input_file.close();
    output_file.close();

    return 0;
}