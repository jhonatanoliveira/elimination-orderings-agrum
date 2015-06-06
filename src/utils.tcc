#ifndef UTILS
#define UTILS

#include <iostream>
#include <vector>

// template<typename int>
void removeInVectorByValue(int element, std::vector<int>& elements) {
    int* p = std::find(elements.begin(), elements.end(), element);
}

#endif