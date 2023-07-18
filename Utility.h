#pragma once
#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <iostream>
#include <vector>
#include <unordered_set>
#include "Struct.h"
#include <tuple>
#include <map>
#include <set>

class Utility {

public:
    Utility();
    ~Utility();
    bool CompareNodeData(const NodeData& a, const NodeData& b);
    std::vector<std::tuple<int, int, int>> removeDuplicates(const std::vector<std::tuple<int, int, int>>& inputVector);
    std::map<std::tuple<int, int, int>, int> findRepeatedTuples(std::vector<std::tuple<int, int, int>>& tuples);
};



#endif // !_UTILITY_
