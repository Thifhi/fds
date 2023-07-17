#include <vector>
#include <cstdint>
#include <cmath>

using namespace std;

// Lookup table data-structure for the upper-bitvector
// Doesn't actually store the upper-bitvector since we only need select queries for Elias-Fano
struct SelectableBitvector {
    vector<uint64_t> table_zeros;
    vector<uint64_t> table_ones;

    SelectableBitvector(vector<bool> data) {
        table_zeros.emplace_back(0);
        table_ones.emplace_back(0);
        for (uint64_t i = 0; i < data.size(); ++i) {
            if (!data[i]) {
                table_zeros.emplace_back(i + 1);
            } else {
                table_ones.emplace_back(i + 1);
            }
        }
    }

    uint64_t select_zero(uint64_t x) {
        if (x >= table_zeros.size()) {
            return UINT64_MAX;
        }
        return table_zeros[x];
    }

    uint64_t select_one(uint64_t x) {
        if (x >= table_zeros.size()) {
            return UINT64_MAX;
        }
        return table_ones[x];
    }

    uint64_t get_size() {
        return (table_zeros.size() + table_ones.size()) * 64;
    }    
};

struct EliasFanoCoded {
    public:
    EliasFanoCoded(vector<uint64_t>& data) : upper_width(set_upper_width(data.size())), lower_width(set_lower_width(data)), upper(build_upper(data)) {
        build_lower(data);
    }

    uint64_t predecessor(uint64_t x) {
        // Extract the lower and upper parts
        uint64_t x_upper = x >> lower_width;
        uint64_t x_lower = x & ((static_cast<uint64_t>(1) << lower_width) - 1);

        // Find the block that contains the possible predecessor in the upper structure
        uint64_t p_candidate = upper.select_zero(x_upper);
        // Find the end of that block so that we can search the indices that lie between for the predecessor
        uint64_t range_high = upper.select_zero(x_upper + 1);
        if (range_high - p_candidate == 1 || range_high == UINT64_MAX) {
            // The block is empty, so return the first number before the block
            if (p_candidate == 0) {
                return UINT64_MAX;
            }
            // Find the first 1 before the query group and get its upper+index, lookup the lower part and sum with the upper
            uint64_t previous_elem_idx = upper.select_one(p_candidate - x_upper);
            uint64_t prev_upper = previous_elem_idx - (p_candidate - x_upper);
            return (prev_upper << lower_width) + lower_at(p_candidate - x_upper - 1);
        }

        // Convert the upper indices to lower indices by subtracting the upper
        p_candidate += -x_upper;
        range_high += -x_upper - 1;

        // Binary search the candidates in the lower
        uint64_t left = p_candidate;
        uint64_t right = range_high - 1;
        uint64_t min_predecessor_lower_idx = UINT64_MAX;

        while (left <= right && right != UINT64_MAX) {
            uint64_t mid = left + (right - left) / 2;

            if (lower_at(mid) <= x_lower) {
                min_predecessor_lower_idx = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }

        if (min_predecessor_lower_idx == UINT64_MAX) {
            // We did not find a smaller or equal value so do the same as we did when the block was empty
            if (p_candidate == 0) {
                return UINT64_MAX;
            }
            uint64_t previous_elem_idx = upper.select_one(p_candidate);
            uint64_t prev_upper = previous_elem_idx - (p_candidate);
            return (prev_upper << lower_width) + lower_at(p_candidate - 1);
        } else {
            return (x_upper << lower_width) + lower_at(min_predecessor_lower_idx);
        }
    }

    uint64_t get_size() {
        return upper.get_size() + lower.size() * 1;
    }

    private:
    uint64_t upper_width;
    uint64_t lower_width;
    SelectableBitvector upper;
    vector<bool> lower;

    static uint64_t set_upper_width(uint64_t n) {
        return ceil(log2(n));
    }

    static uint64_t set_lower_width(vector<uint64_t>& data) {
        // Set the universe size to the largest number in data + 1 to use as little bits as possible
        uint64_t max = 0;
        for (auto item : data) {
            if (item > max) {
                max = item;
            }
        }
        return ceil(log2(max + 1) - log2(data.size()));
    }

    // Builds the upper-bitvector of Elias-Fano
    vector<bool> build_upper(vector<uint64_t>& data) {
        vector<bool> bitvector(data.size() + pow(2, upper_width) + 1);
        for (uint64_t i = 0; i < data.size(); ++i) {
            uint64_t bv_idx = (data[i] >> lower_width) + i;
            bitvector[bv_idx] = 1;
        }
        return bitvector;
    }

    // Builds the lower-bitvector of Elias-Fano
    void build_lower(vector<uint64_t>& data) {
        lower.resize(data.size() * lower_width);
        for (uint64_t i = 0; i < data.size(); ++i) {
            for (uint64_t j = 0; j < lower_width; ++j) {
                lower[i * lower_width + j] = (data[i] >> (lower_width - j - 1)) & 1;
            }
        }
    }

    // Returns the lower-part number at index n by converting the bits in the bitvector into a UINT64_T
    uint64_t lower_at(uint64_t n) {
        uint64_t res = 0;
        for (uint64_t i = 0; i < lower_width; ++i) {
            res = res << 1;
            res += lower[lower_width * n + i];
        }
        return res;
    }
};