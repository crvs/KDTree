/*
 * file: KDTree.hpp
 * author: J. Frederico Carvalho
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 * https://rosettacode.org/wiki/K-d_tree
 *
 * It is a reimplementation of the C code using C++.  It also includes a few
 * more queries than the original, namely finding all points at a distance
 * smaller than some given distance to a point.
 *
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>

#include "KDTree.hpp"

KDNode::KDNode() = default;

KDNode::KDNode(point_t const& pt, size_t const& idx_, KDNodePtr const& left_,
               KDNodePtr const& right_) {
    x = pt;
    index = idx_;
    left = left_;
    right = right_;
}

KDNode::KDNode(pointIndex const& pi, KDNodePtr const& left_,
               KDNodePtr const& right_) {
    x = pi.first;
    index = pi.second;
    left = left_;
    right = right_;
}

KDNode::~KDNode() = default;

double KDNode::coord(size_t const& idx) { return x.at(idx); }
KDNode::operator bool() { return (!x.empty()); }
KDNode::operator point_t() { return x; }
KDNode::operator size_t() { return index; }
KDNode::operator pointIndex() { return pointIndex(x, index); }

KDNodePtr NewKDNodePtr() {
    KDNodePtr mynode = std::make_shared<KDNode>();
    return mynode;
}

inline double dist2(point_t const& a, point_t const& b) {
    double distc = 0;
    for (size_t i = 0; i < a.size(); i++) {
        double di = a.at(i) - b.at(i);
        distc += di * di;
    }
    return distc;
}

inline double dist2(KDNodePtr const& a, KDNodePtr const& b) {
    return dist2(a->x, b->x);
}

inline double dist(point_t const& a, point_t const& b) {
    return std::sqrt(dist2(a, b));
}

inline double dist(KDNodePtr const& a, KDNodePtr const& b) {
    return std::sqrt(dist2(a, b));
}

comparer::comparer(size_t idx_) : idx{idx_} {};

inline bool comparer::compare_idx(pointIndex const& a, pointIndex const& b) {
    return (a.first.at(idx) < b.first.at(idx));
}

inline void sort_on_idx(pointIndexArr::iterator const& begin,
                        pointIndexArr::iterator const& end, size_t idx) {
    comparer comp(idx);
    comp.idx = idx;

    using std::placeholders::_1;
    using std::placeholders::_2;

    std::nth_element(begin, begin + std::distance(begin, end) / 2, end,
                     std::bind(&comparer::compare_idx, comp, _1, _2));
}

using pointVec = std::vector<point_t>;

KDNodePtr KDTree::make_tree(pointIndexArr::iterator const& begin,
                            pointIndexArr::iterator const& end,
                            size_t const& length, size_t const& level) {
    if (begin == end) {
        return NewKDNodePtr(); // empty tree
    }

    size_t dim = begin->first.size();

    if (length > 1) {
        sort_on_idx(begin, end, level);
    }

    auto middle = begin + (length / 2);

    auto l_begin = begin;
    auto l_end = middle;
    auto r_begin = middle + 1;
    auto r_end = end;

    size_t l_len = length / 2;
    size_t r_len = length - l_len - 1;

    KDNodePtr left;
    if (l_len > 0 && dim > 0) {
        left = make_tree(l_begin, l_end, l_len, (level + 1) % dim);
    } else {
        left = leaf_;
    }
    KDNodePtr right;
    if (r_len > 0 && dim > 0) {
        right = make_tree(r_begin, r_end, r_len, (level + 1) % dim);
    } else {
        right = leaf_;
    }

    // KDNode result = KDNode();
    return std::make_shared<KDNode>(*middle, left, right);
}

KDTree::KDTree(pointVec point_array) : leaf_{std::make_shared<KDNode>()} {
    // iterators
    pointIndexArr arr;
    for (size_t i = 0; i < point_array.size(); i++) {
        arr.push_back(pointIndex(point_array.at(i), i));
    }

    auto begin = arr.begin();
    auto end = arr.end();

    size_t length = arr.size();
    size_t level = 0; // starting

    root_ = KDTree::make_tree(begin, end, length, level);
}

KDNodePtr KDTree::nearest_(KDNodePtr const& branch, point_t const& pt,
                           size_t const& level, KDNodePtr const& best,
                           double const& best_dist) {
    double d, dx, dx2;

    if (!bool(*branch)) {
        return NewKDNodePtr(); // basically, null
    }

    point_t branch_pt(*branch);
    size_t dim = branch_pt.size();

    d = dist2(branch_pt, pt);
    dx = branch_pt.at(level) - pt.at(level);
    dx2 = dx * dx;

    KDNodePtr best_l = best;
    double best_dist_l = best_dist;

    if (d < best_dist) {
        best_dist_l = d;
        best_l = branch;
    }

    size_t next_lv = (level + 1) % dim;
    KDNodePtr section;
    KDNodePtr other;

    // select which branch makes sense to check
    if (dx > 0) {
        section = branch->left;
        other = branch->right;
    } else {
        section = branch->right;
        other = branch->left;
    }

    // keep nearest neighbor from further down the tree
    KDNodePtr further = nearest_(section, pt, next_lv, best_l, best_dist_l);
    if (!further->x.empty()) {
        double dl = dist2(further->x, pt);
        if (dl < best_dist_l) {
            best_dist_l = dl;
            best_l = further;
        }
    }
    // only check the other branch if it makes sense to do so
    if (dx2 < best_dist_l) {
        further = nearest_(other, pt, next_lv, best_l, best_dist_l);
        if (!further->x.empty()) {
            double dl = dist2(further->x, pt);
            if (dl < best_dist_l) {
                best_dist_l = dl;
                best_l = further;
            }
        }
    }

    return best_l;
};

// default caller
KDNodePtr KDTree::nearest_(point_t const& pt) {
    size_t level = 0;
    double branch_dist = dist2(point_t(*root_), pt);
    return nearest_(root_,        // beginning of tree
                    pt,           // point we are querying
                    level,        // start from level 0
                    root_,        // best is the root
                    branch_dist); // best_dist = branch_dist
};

point_t KDTree::nearest_point(point_t const& pt) {
    return point_t(*nearest_(pt));
};

size_t KDTree::nearest_index(point_t const& pt) {
    return size_t(*nearest_(pt));
};

pointIndex KDTree::nearest_pointIndex(point_t const& pt) {
    KDNodePtr Nearest = nearest_(pt);
    return pointIndex(point_t(*Nearest), size_t(*Nearest));
}

void KDTree::neighborhood_(KDNodePtr const& branch, point_t const& pt,
                           double const& rad, size_t const& level,
                           pointIndexArr& nbh) {
    double d, dx, dx2;

    if (!bool(*branch)) {
        // branch has no point, means it is a leaf,
        // no points to add
        return;
    }

    size_t dim = pt.size();

    double r2 = rad * rad;

    d = dist2(point_t(*branch), pt);
    dx = point_t(*branch).at(level) - pt.at(level);
    dx2 = dx * dx;

    if (d <= r2) {
        nbh.push_back(pointIndex(*branch));
    }

    KDNodePtr section;
    KDNodePtr other;
    if (dx > 0) {
        section = branch->left;
        other = branch->right;
    } else {
        section = branch->right;
        other = branch->left;
    }

    neighborhood_(section, pt, rad, (level + 1) % dim, nbh);
    if (dx2 < r2) {
        neighborhood_(other, pt, rad, (level + 1) % dim, nbh);
    }
};

pointIndexArr KDTree::neighborhood(point_t const& pt, double const& rad) {
    size_t level = 0;
    pointIndexArr nbh;
    neighborhood_(root_, pt, rad, level, nbh);
    return nbh;
}

pointVec KDTree::neighborhood_points(point_t const& pt, double const& rad) {
    size_t level = 0;
    auto nbh = std::make_shared<pointIndexArr>();
    neighborhood_(root_, pt, rad, level, *nbh);
    pointVec nbhp;
    nbhp.resize(nbh->size());
    std::transform(nbh->begin(), nbh->end(), nbhp.begin(),
                   [](pointIndex x) { return x.first; });
    return nbhp;
}

indexArr KDTree::neighborhood_indices(point_t const& pt, double const& rad) {
    size_t level = 0;
    auto nbh = std::make_shared<pointIndexArr>();
    neighborhood_(root_, pt, rad, level, *nbh);
    indexArr nbhi;
    nbhi.resize(nbh->size());
    std::transform(nbh->begin(), nbh->end(), nbhi.begin(),
                   [](pointIndex x) { return x.second; });
    return nbhi;
}
