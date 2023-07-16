#include <vector>
#include <cstdint>
#include <cmath>

using namespace std;

struct SelectableBitvector {
    vector<bool> data;
    vector<uint64_t> table_zeros;
    vector<uint64_t> table_ones;

    SelectableBitvector(vector<bool> data) : data(data) {
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
};

struct EliasFanoCoded {
    public:
    EliasFanoCoded(vector<uint64_t>& data) : upper_width(set_upper_width(data.size())), lower_width(set_lower_width(data)), upper(build_upper(data)) {
        build_lower(data);
    }

    uint64_t predecessor_idx(uint64_t x) {
        uint64_t x_upper = x >> lower_width;
        uint64_t x_lower = x & ((static_cast<uint64_t>(1) << lower_width) - 1);
        uint64_t p_candidate = upper.select_zero(x_upper);
        uint64_t range_high = upper.select_zero(x_upper + 1);
        if (range_high - p_candidate == 1 || range_high == UINT64_MAX) {
            if (p_candidate == 0) {
                return UINT64_MAX;
            }
            uint64_t previous_elem_idx = upper.select_one(p_candidate - x_upper);
            uint64_t prev_upper = previous_elem_idx - (p_candidate - x_upper);
            return (prev_upper << lower_width) + lower_at(p_candidate - x_upper - 1);
        }

        p_candidate += -x_upper;
        range_high += -x_upper - 1;

        // Binary search the candidates
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
            // We did not find a smaller or equal value
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

    private:
    uint64_t upper_width;
    uint64_t lower_width;
    SelectableBitvector upper;
    vector<bool> lower;

    static uint64_t set_upper_width(uint64_t n) {
        return ceil(log2(n));
    }

    static uint64_t set_lower_width(vector<uint64_t>& data) {
        uint64_t max = 0;
        for (auto item : data) {
            if (item > max) {
                max = item;
            }
        }
        return ceil(log2(max + 1) - log2(data.size()));
    }

    vector<bool> build_upper(vector<uint64_t>& data) {
        vector<bool> bitvector(data.size() + pow(2, upper_width) + 1);
        for (uint64_t i = 0; i < data.size(); ++i) {
            uint64_t bv_idx = (data[i] >> lower_width) + i;
            bitvector[bv_idx] = 1;
        }
        return bitvector;
    }

    void build_lower(vector<uint64_t>& data) {
        lower.resize(data.size() * lower_width);
        for (uint64_t i = 0; i < data.size(); ++i) {
            for (uint64_t j = 0; j < lower_width; ++j) {
                lower[i * lower_width + j] = (data[i] >> (lower_width - j - 1)) & 1;
            }
        }
    }

    uint64_t lower_at(uint64_t n) {
        uint64_t res = 0;
        for (uint64_t i = 0; i < lower_width; ++i) {
            res = res << 1;
            res += lower[lower_width * n + i];
        }
        return res;
    }
};