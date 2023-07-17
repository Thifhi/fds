#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>

#include "rmq/api.h"
#include "predecessor/api.h"

using namespace std;

vector<uint64_t> read_data(ifstream& input_file) {
    int n;
    input_file >> n;
    vector<uint64_t> data(n);

    for (int i = 0; i < data.size(); ++i) {
        uint64_t num;
        input_file >> num;
        data[i] = num;
    }
    
    return data;
}

vector<pair<uint64_t, uint64_t>> read_rmq_queries(ifstream& input_file) {
    vector<pair<uint64_t, uint64_t>> queries;

    char comma;
    uint64_t s;
    uint64_t e;
    while (input_file >> s) {
        input_file >> comma;
        input_file >> e;
        queries.emplace_back(make_pair(s, e));
    }
    
    return queries;
}

vector<uint64_t> read_pd_queries(ifstream& input_file) {
    vector<uint64_t> queries;

    uint64_t query;
    while (input_file >> query) {
        queries.emplace_back(query);
    }

    return queries;
}

void output_results(tuple<vector<uint64_t>, int64_t, uint64_t>& results, ofstream& output_file, string algorithm) {
    for (auto result : get<0>(results)) {
        output_file << result << endl;
    }
    output_file << "RESULT algo=" << algorithm << " name=mustafa_enes_batur time=" << get<1>(results) << " space=" << get<2>(results);
}

int main(int argc, char* argv[]) {
    string algorithm = argv[1];
    string input_file_name = argv[2];

    // Open the input file
    ifstream input_file(input_file_name);
    if (!input_file) {
        std::cerr << "Failed to open input file: " << input_file_name << endl;
        return 1;
    }

    string output_file_name = argv[3];
    // Open the output file
    ofstream output_file(output_file_name);
    if (!output_file) {
        cerr << "Failed to open output file: " << output_file_name << endl;
        input_file.close();
        return 1;
    }

    auto data = read_data(input_file);
    tuple<vector<uint64_t>, int64_t, uint64_t> results;

    if (algorithm == "pd") {
        auto queries = read_pd_queries(input_file);
        results = solve_pd(data, queries);
    } else if (algorithm == "rmq") {
        auto queries = read_rmq_queries(input_file);
        results = solve_rmq(data, queries, "linear");
    } else {
        cout << "Unknown algorithm!" << endl;
        exit(1);
    }

    output_results(results, output_file, algorithm);

    // Close the files
    input_file.close();
    output_file.close();

    return 0;
}