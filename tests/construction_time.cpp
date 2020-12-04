#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <ctime>
#include <ratio>
#include <numeric>
#include <algorithm>
#include "KDTree.hpp"

#define DIM 3

double getNum()
{
	return ((double)rand() / (RAND_MAX));
}

std::vector<double> generateVector()
{
    std::vector<double> temp(DIM);
	for (size_t idx = 0; idx < DIM; idx++)
    {
		temp[idx] = getNum();
    }
	return temp;
}

std::vector<std::vector<double>> getListofGeneratedVectors(size_t length)
{
    std::vector<std::vector<double>> temp(length);
	for (size_t idx = 0; idx < length; idx++)
    {
		temp[idx] = generateVector();
    }
	return temp;
}

snt main()
{
    // seed
	srand(5);

    size_t npoints = 400000;
    std::cout << "constructing KDTree with " << npoints << " points." << std::endl;

    td::vector<point_t> points = getListofGeneratedVectors(npoints);

    uto start = std::chrono::high_resolution_clock::now();
    KDTree tree(points);
    auto stop = std::chrono::high_resolution_clock::now();
    auto timespan = std::chrono::duration<double>(stop - start);
    std::cout << "it took " << timespan.count() << " seconds.";
    return 0;
}
