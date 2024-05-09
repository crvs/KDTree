#include "KDTree.hpp"
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>
#include <numeric>
#include <random>
#include <ratio>
#include <vector>

#define DIM 3

double getNum() { return ((double)rand() / (RAND_MAX)); }

std::vector<double> generateVector() {
    std::vector<double> temp(DIM);
    for (size_t idx = 0; idx < DIM; idx++) {
        temp[idx] = getNum();
    }
    return temp;
}

std::vector<std::vector<double>> getListofGeneratedVectors(size_t length) {
    std::vector<std::vector<double>> temp(length);
    for (size_t idx = 0; idx < length; idx++) {
        temp[idx] = generateVector();
    }
    return temp;
}

int main() {
    // seed
    srand(5);

    size_t npoints = 4'000'000;
    std::cout << "constructing KDTree with " << npoints << " points..."
              << std::endl;

    std::vector<point_t> points = getListofGeneratedVectors(npoints);

    auto start = std::chrono::high_resolution_clock::now();
    KDTree tree(points);
    auto stop = std::chrono::high_resolution_clock::now();
    auto timespan = std::chrono::duration<double>(stop - start);
    std::cout << "It took " << timespan.count() << " seconds." << std::endl;
    return 0;
}
