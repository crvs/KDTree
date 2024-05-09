// @Authors: Sushil B. and Paul M. - (C) 2019, MIT License

#include "KDTree.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#define DIM 3

double getNum() { return ((double)rand() / (RAND_MAX)); }
using secondsf = std::chrono::duration<float>;

std::vector<double> generateVector() {
    std::vector<double> temp(DIM);
    for (size_t idx = 0; idx < DIM; idx++) {
        temp[idx] = getNum();
    }
    return temp;
}

std::vector<std::vector<double>> getListofGeneratedVectors(int length) {
    std::vector<std::vector<double>> temp(length);
    for (size_t idx = 0; idx < static_cast<size_t>(length); idx++) {
        temp[idx] = generateVector();
    }
    return temp;
}

double sumSqrdErr(std::vector<double>& p1, std::vector<double>& p2) {
    std::vector<double> diff(DIM);
    std::vector<double> square(DIM);
    std::transform(p1.begin(), p1.end(), p2.begin(), diff.begin(),
                   std::minus<double>());
    std::transform(diff.begin(), diff.end(), diff.begin(), square.begin(),
                   std::multiplies<double>());
    return std::accumulate(square.begin(), square.end(), 0.0);
}

int main() {
    // seed
    srand(5);

    const std::vector<int> dataPointSizes = {50,  100, 200, 300, 400,  500,
                                             600, 700, 800, 900, 1000, 2000};
    constexpr int nIter = 50;

    std::cout << "Total number of iterations ran: " << nIter << std::endl;

    for (auto& sizes : dataPointSizes) {
        std::vector<std::vector<double>> points(sizes,
                                                std::vector<double>(DIM)),
            pointToRetrieve(sizes, std::vector<double>(DIM));
        std::chrono::nanoseconds kdTreeRetTotalTime{};
        std::chrono::nanoseconds bruteForceRetTotalTime{};

        int correct = 0;
        for (int i = 0; i < nIter; i++) {
            // generate test points to build a tree
            points = getListofGeneratedVectors(sizes);

            // genereate KD Tree
            KDTree tree(points);

            // generate retrieve test data points
            pointToRetrieve = getListofGeneratedVectors(sizes);

            for (auto& vals : pointToRetrieve) {
                double minSumSqdErr = std::numeric_limits<double>::max();
                std::vector<double> groundTruthVec(DIM);

                auto bf_start = std::chrono::high_resolution_clock::now();
                for (auto& gtvals : points) {
                    double sumSqdErr = sumSqrdErr(gtvals, vals);
                    if (sumSqdErr < minSumSqdErr) {
                        minSumSqdErr = sumSqdErr;
                        groundTruthVec = gtvals;
                    }
                }
                bruteForceRetTotalTime +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() - bf_start);

                std::vector<double> checkVec(DIM);

                auto kdt_start = std::chrono::high_resolution_clock::now();
                checkVec = tree.nearest_point(vals);
                kdTreeRetTotalTime +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() - kdt_start);

                if (std::equal(groundTruthVec.begin(), groundTruthVec.end(),
                               checkVec.begin())) {
                    correct += 1;
                }
            }
        }
        std::cout << "Accuracy (tested with " << sizes
                  << " datasets per iter): "
                  << ((correct * 100.0) / (sizes * nIter))
                  << "%. Total Number of correct queries: "
                  << (int)(correct / nIter) << " / " << sizes << std::endl;
        if (correct != nIter * sizes) {
            return 1;
        }
        std::cout
            << "Total query time: { bruteForce: "
            << std::chrono::duration_cast<secondsf>(bruteForceRetTotalTime)
                   .count()
            << ", kdTree: "
            << std::chrono::duration_cast<secondsf>(kdTreeRetTotalTime).count()
            << " }" << std::endl;
    }
    return 0;
}
