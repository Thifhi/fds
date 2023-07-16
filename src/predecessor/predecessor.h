#include <vector>
#include <cstdint>
#include <cmath>

using namespace std;

struct SelectableBitvector {
    vector<bool> data;
    vector<uint64_t> table;

    SelectableBitvector(vector<bool> data) : data(data) {
        table.resize(data.size());

        uint64_t found = 0;
        table[0] = 0;
        for (uint64_t i = 0; i < data.size(); ++i) {
            if (!data[i]) {
                ++found;
                table[found] = i + 1;
            }
        }
    }

    uint64_t select_zero(uint64_t x) {
        return table[x];
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
        if (p_candidate - range_high == 1 || range_high == UINT64_MAX) {
            // No elements in the block
            return (p_candidate - x_upper - 1);
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
            return p_candidate - 1;
        } else {
            return min_predecessor_lower_idx;
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