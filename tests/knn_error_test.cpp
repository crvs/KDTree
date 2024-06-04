// @Authors: Sushil B. and Paul M. - (C) 2019, MIT License

#include "KDTree.hpp"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#define DIM 2

double getNum() { return ((double)rand() / (RAND_MAX)); }
using secondsf = std::chrono::duration<float>;
using pointDistance = std::pair<std::vector<double>, double>;
bool samePt(std::vector<double> v1, std::vector<double> v2) {
    auto const result = std::equal(v1.begin(), v1.end(), v2.begin());
    if (!result && false) {
        std::cout << "Different points: [ ";
        for (auto const v : v1) {
            std::cout << v << " ";
        }
        std::cout << "] [ ";
        for (auto const v : v2) {
            std::cout << v << " ";
        }
        std::cout << "]" << std::endl;
    }
    return result;
}

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

double sumSqrdErr(std::vector<double> const& p1,
                  std::vector<double> const& p2) {
    std::vector<double> diff(DIM);
    std::vector<double> square(DIM);
    std::transform(p1.begin(), p1.end(), p2.begin(), diff.begin(),
                   std::minus<double>());
    std::transform(diff.begin(), diff.end(), diff.begin(), square.begin(),
                   std::multiplies<double>());
    return std::accumulate(square.begin(), square.end(), 0.0);
}

int main() {
    bool success = true;
    // seed
    srand(5);

    const std::vector<int> dataPointSizes = {50, 500, 1000};
    constexpr int nIter = 30;

    std::cout << "Total number of iterations ran: " << nIter << std::endl;

    for (auto& sizes : dataPointSizes) {
        std::vector<std::vector<double>> points(sizes,
                                                std::vector<double>(DIM)),
            pointToRetrieve(sizes, std::vector<double>(DIM));
        std::chrono::nanoseconds kdTreeRetTotalTime{};
        std::chrono::nanoseconds bruteForceRetTotalTime{};

        int correct = 0;
        int const k = 10;
        for (int i = 0; i < nIter; i++) {
            // generate test points to build a tree
            points = getListofGeneratedVectors(sizes);

            // genereate KD Tree
            KDTree tree(points);

            // generate retrieve test data points
            pointToRetrieve = getListofGeneratedVectors(sizes);

            for (auto& vals : pointToRetrieve) {
                std::vector<double> groundTruthVec(DIM);

                auto bf_start = std::chrono::high_resolution_clock::now();
                std::list<pointDistance> point_distances{};
                for (auto& gtvals : points) {
                    pointDistance const point_distance =
                        std::make_pair(gtvals, sumSqrdErr(gtvals, vals));
                    point_distances.insert(
                        std::upper_bound(point_distances.begin(),
                                         point_distances.end(), point_distance,
                                         [](auto const& a, auto const& b) {
                                             return a.second < b.second;
                                         }),
                        point_distance);
                    if (point_distances.size() > k) {
                        point_distances.pop_back();
                    }
                }
                bruteForceRetTotalTime +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() - bf_start);

                auto kdt_start = std::chrono::high_resolution_clock::now();
                auto const checkVec = tree.nearest_points(vals, k);
                kdTreeRetTotalTime +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() - kdt_start);

                if (checkVec.size() == k &&
                    std::equal(point_distances.begin(), point_distances.end(),
                               checkVec.begin(),
                               [](auto const& v1, auto const& v2) {
                                   return samePt(v1.first, v2);
                               })) {
                    correct += 1;
                } else {
                    // assert(false);
                    std::cout << "Distance vectors: gt: [ ";
                    for (auto const& pd : point_distances) {
                        std::cout << std::log(pd.second) << " ";
                    }
                    std::cout << "] pred: [ ";
                    for (auto const& v : checkVec) {
                        std::cout << std::log(sumSqrdErr(v, vals)) << " ";
                    }
                    std::cout << "]" << std::endl;
                    success = false;
                }
            }
        }
        std::cout << "Accuracy (tested with " << sizes
                  << " datasets per iter): "
                  << ((correct * 100.0) / (sizes * nIter))
                  << "%. Total Number of correct queries: "
                  << (int)(correct / nIter) << " / " << sizes << std::endl;
        std::cout
            << "Total query time: { bruteForce: "
            << std::chrono::duration_cast<secondsf>(bruteForceRetTotalTime)
                   .count()
            << ", kdTree: "
            << std::chrono::duration_cast<secondsf>(kdTreeRetTotalTime).count()
            << " }" << std::endl;
    }
    if (success)
        return 0;
    return 1;
}
