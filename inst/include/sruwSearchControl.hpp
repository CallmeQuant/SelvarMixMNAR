#ifndef SRUW_SEARCH_CONTROL_H
#define SRUW_SEARCH_CONTROL_H

#include <algorithm>
#include <cstddef>
#include <vector>

inline int sruw_regression_iteration_cap(std::size_t predictorCount) {
    const long long p = static_cast<long long>(predictorCount);
    const long long scaled = 4LL * p * p + 16LL;
    return static_cast<int>(std::min(10000LL, std::max(16LL, scaled)));
}

inline std::vector<int> sruw_regression_state_key(
    const std::vector<int>& selected,
    const std::vector<int>& lastIncluded,
    const std::vector<int>& lastExcluded) {
    std::vector<int> key = selected;
    std::sort(key.begin(), key.end());
    key.push_back(0);

    std::vector<int> included = lastIncluded;
    std::sort(included.begin(), included.end());
    key.insert(key.end(), included.begin(), included.end());
    key.push_back(0);

    std::vector<int> excluded = lastExcluded;
    std::sort(excluded.begin(), excluded.end());
    key.insert(key.end(), excluded.begin(), excluded.end());
    return key;
}

#endif
