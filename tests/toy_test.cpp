#include <iostream>
#include <vector>

#include "KDTree.hpp"

// using point_t = std::vector< double >;
// using pointVec = std::vector< point_t >;

// point_t pt(2);

int main() {
    pointVec points;
    point_t pt;

    pt = {0.0, 0.0};
    points.push_back(pt);
    pt = {1.0, 0.0};
    points.push_back(pt);
    pt = {0.0, 1.0};
    points.push_back(pt);
    pt = {1.0, 1.0};
    points.push_back(pt);
    pt = {0.5, 0.5};
    points.push_back(pt);

    KDTree tree(points);

    std::cout << "nearest test\n";
    pt = {0.8, 0.2};
    auto res = tree.nearest_point(pt);
    for (double b : res) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    /*
    std::cout << "going down the tree\n";

    for (auto b : point_t(*tree)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    for (auto b : point_t(*tree->left)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    for (auto b : point_t(*tree->right)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    for (auto b : point_t(*tree->left->left)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    for (auto b : point_t(*tree->right->left)) {
        std::cout << b << " ";
    }
    std::cout << '\n';

    std::cout << "printing nbh\n";

    pt = {.0, .5};

    */
    auto res2 = tree.neighborhood_points(pt, .55);

    for (point_t a : res2) {
        for (double b : a) {
            std::cout << b << " ";
        }
        std::cout << '\n';
    }
    return 0;
}
