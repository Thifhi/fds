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

void output_results(vector<uint64_t>& results, ofstream& output_file) {

}

int main(int argc, char* argv[]) {
    string algorithm = argv[1];
    string inputFileName = argv[2];

    // Open the input file
    ifstream input_file(inputFileName);
    if (!input_file) {
        std::cerr << "Failed to open input file: " << inputFileName << endl;
        return 1;
    }

    string outputFileName = argv[3];
    // Open the output file
    ofstream output_file(outputFileName);
    if (!output_file) {
        cerr << "Failed to open output file: " << outputFileName << endl;
        input_file.close();
        return 1;
    }

    auto data = read_data(input_file);
    vector<uint64_t> results;

    if (algorithm == "pd") {
        auto queries = read_pd_queries(input_file);
        results = solve_pd(data, queries);
    } else if (algorithm == "rmq") {
        auto queries = read_rmq_queries(input_file);
        results = solve_rmq(data, queries);
    } else {
        cout << "Unknown algorithm!" << endl;
        exit(1);
    }

    output_results(results, output_file);

    // Close the files
    input_file.close();
    output_file.close();

    return 0;
}