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
using pointDistance = std::pair<std::vector<double>, double>;
bool samePt(std::vector<double> v1, std::vector<double> v2) {
    auto const result = std::equal(v1.begin(), v1.end(), v2.begin(),
                                   [](auto const& p1, auto const& p2) {
                                       std::cout << std::abs(p1 - p2) << " ";
                                       return std::abs(p1 - p2) < 0.1;
                                   });
    std::cout << std::endl;
    if (!result) {
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
        int const k = 4;
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
                    double sumSqdErr = sumSqrdErr(gtvals, vals);
                    pointDistance const point_distance =
                        std::make_pair(gtvals, sumSqdErr);
                    point_distances.insert(
                        std::upper_bound(point_distances.begin(),
                                         point_distances.end(), point_distance,
                                         [](auto const& a, auto const& b) {
                                             return a.second < b.second;
                                         }),
                        point_distance);
                    if (point_distances.size() > 10) {
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
