#include <vector>
#include <cstdint>
#include <cmath>

using namespace std;

struct SelectableBitvector {
    vector<bool> data;

    SelectableBitvector(vector<bool> data) : data(data) {}

    uint64_t select_zero(uint64_t x) {
        if (x == 0) {
            return 0;
        }

        uint64_t found = 0;
        for (uint64_t i = 0; i < data.size(); ++i) {
            if (!data[i]) {
                ++found;
            }
            if (found == x) {
                return i + 1;
            }
        }
        return UINT64_MAX;
    }
};

struct EliasFanoCoded {
    public:
    EliasFanoCoded(vector<uint64_t>& data) : upper_width(set_upper_width(data.size())), lower_width(set_lower_width(data.size())), upper(build_upper(data)) {
        build_lower(data);
    }

    uint64_t predecessor_idx(uint64_t x) {
        uint64_t x_upper = x >> lower_width;
        uint64_t x_lower = x & ((1 << lower_width) - 1);
        uint64_t p_candidate = upper.select_zero(x_upper);
        uint64_t range_high = upper.select_zero(x_upper + 1);
        if (p_candidate - range_high == 1 || range_high == UINT64_MAX) {
            // No elements in the block
            return max(p_candidate - x_upper - 1, static_cast<uint64_t>(0));
        }

        p_candidate += -x_upper;
        range_high += -x_upper - 1;

        // Binary search the candidates
        uint64_t left = p_candidate;
        uint64_t right = range_high - 1;
        uint64_t min_predecessor_lower_idx = UINT64_MAX;

        while (left <= right) {
            int mid = left + (right - left) / 2;

            if (lower_at(mid) <= x_lower) {
                min_predecessor_lower_idx = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }

        if (min_predecessor_lower_idx == UINT64_MAX) {
            // We did not find a smaller or equal value
            return max(p_candidate - 1, static_cast<uint64_t>(0));
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

    static uint64_t set_lower_width(uint64_t n) {
        return 2;
        // return ceil(log2(UINT64_MAX) - log2(n));
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
                lower[i * lower_width + j] = (data[i] >> j) & 1;
            }
        }
    }

    uint64_t lower_at(uint64_t n) {
        uint64_t res = 0;
        for (uint64_t i = 0; i < lower_width; ++i) {
            res = res << 1; // First iteration doesn't matter since we are at 0
            res = res | (lower[lower_width * n + i]);
        }
        return res;
    }
};