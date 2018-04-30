#pragma once
/*
 * file: KDTree.hpp
 * author: J. Frederico Carvalho
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 *  https://rosettacode.org/wiki/K-d_tree
 * It is a reimplementation of the C code using C++.
 * It also includes a few more queries than the original
 *
 */

/*
 * TODO:
 * - create a KDTree class and wrap the basic queries in it
 * - make output of basic queries in terms of point_t rather than KDNodePtr
 */

#include <algorithm>
#include <functional>

#include <iterator>
#include <memory>
#include <vector>

using point_t = std::vector< double >;

class KDNode {
   public:
    using KDNodePtr = std::shared_ptr< KDNode >;

    point_t x;
    KDNodePtr left, right;

    KDNode() { x = point_t(); }

    KDNode(const point_t &pt, const KDNodePtr &left_, const KDNodePtr &right_) {
        x = pt;
        left = left_;
        right = right_;
    }

    double coord(size_t idx) { return x.at(idx); }

    // conversions
    explicit operator bool() { return (!x.empty()); };
    explicit operator point_t() { return x; };
};

using KDNodePtr = std::shared_ptr< KDNode >;

KDNodePtr NewKDNodePtr() { return std::make_shared< KDNode >(); }

using KDNodePArray = std::vector< KDNodePtr >;

inline double dist(const point_t &a, const point_t &b) {
    double distc = 0;
    for (size_t i = 0; i < a.size(); i++) {
        double di = a.at(i) - b.at(i);
        distc += di * di;
    }
    return distc;
}

inline double dist(const KDNodePtr &a, const KDNodePtr &b) {
    return dist(a->x, b->x);
}

template < typename T >
class comparer {
   public:
    size_t idx;

    comparer() = default;

    explicit comparer(size_t idx_) : idx{idx_} {};

    inline bool compare_idx(const std::vector< T > &a,  //
                            const std::vector< T > &b   //
    ) {
        return (a.at(idx) < b.at(idx));  //
    }
};

/* see quickselect method */
template < typename T >
inline void sort_on_idx(
    typename std::vector< std::vector< T > >::iterator &begin,
    typename std::vector< std::vector< T > >::iterator &end,  //
    size_t idx) {
    comparer< double > c(idx);
    c.idx = idx;

    using std::placeholders::_1;
    using std::placeholders::_2;

    std::sort(begin, end,
              std::bind(&comparer< double >::compare_idx, c, _1, _2));
}

using pointVec = std::vector< point_t >;

class KDTree {
    KDNodePtr root;

    static KDNodePtr make_tree(pointVec::iterator begin,  //
                               pointVec::iterator end,    //
                               size_t length,             //
                               size_t level               //
    ) {
        if (begin == end) {
            return NewKDNodePtr();  // empty tree
        }

        size_t dim = begin->size();

        sort_on_idx< double >(begin, end, level);
        auto middle = begin + (length / 2);

        auto l_begin = begin;
        auto l_end = middle;
        auto r_begin = middle + 1;
        auto r_end = end;

        size_t l_len = length / 2;
        size_t r_len = (length / 2) + (length % 2);

        KDNodePtr left = make_tree(l_begin, l_end, l_len, (level + 1) % dim);
        KDNodePtr right = make_tree(r_begin, r_end, r_len, (level + 1) % dim);

        // KDNode result = KDNode();
        return std::make_shared< KDNode >(*middle, left, right);
    }

   public:
    explicit KDTree(pointVec point_array) {
        // iterators
        auto begin = point_array.begin();
        auto end = point_array.end();

        size_t length = point_array.size();
        size_t level = 0;  // starting

        root = KDTree::make_tree(begin, end, length, level);
    }

   private:
    KDNodePtr nearest_(           //
        const KDNodePtr &branch,  //
        const point_t &pt,        //
        const size_t &level,      //
        const KDNodePtr &best,    //
        const double &best_dist   //
    ) {
        double d, dx, dx2;

        if (!bool(*branch)) {
            return NewKDNodePtr();  // basically, null
        }

        point_t branch_pt(*branch);
        size_t dim = branch_pt.size();

        d = dist(branch_pt, pt);
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
            double dl = dist(further->x, pt);
            if (dl < best_dist_l) {
                best_dist_l = dl;
                best_l = further;
            }
            // only check the other branch if it makes sense to do so
            if (dx2 < best_dist_l) {
                further = nearest_(other, pt, next_lv, best_l, best_dist_l);
                if (!further->x.empty()) {
                    dl = dist(further->x, pt);
                    if (dl < best_dist_l) {
                        best_dist_l = dl;
                        best_l = further;
                    }
                }
            }
        }

        return best_l;
    };

    // default caller
    KDNodePtr nearest_(const KDNodePtr &branch, const point_t &pt) {
        size_t level = 0;
        // KDNodePtr best = branch;
        double branch_dist = dist(point_t(*branch), pt);
        return nearest_(branch,        // beginning of tree
                        pt,            // point we are querying
                        level,         // start from level 0
                        branch,        // best = branch
                        branch_dist);  // best_dist = branch_dist
    }

   public:
    point_t nearest(const point_t &pt) { return point_t(*nearest_(root, pt)); }

   private:
    pointVec neighborhood_(       //
        const KDNodePtr &branch,  //
        const point_t &pt,        //
        const double &rad,        //
        const size_t &level       //
    ) {
        double d, dx, dx2;

        if (!bool(*branch)) {
            // branch has no point, means it is a leaf,
            // no points to add
            return pointVec();
        }

        size_t dim = pt.size();

        double r2 = rad * rad;

        d = dist(point_t(*branch), pt);
        dx = point_t(*branch).at(level) - pt.at(level);
        dx2 = dx * dx;

        pointVec nbh, nbh_s, nbh_o;
        if (d <= r2) {
            nbh.push_back(point_t(*branch));
        }

        //
        KDNodePtr section;
        KDNodePtr other;
        if (dx > 0) {
            section = branch->left;
            other = branch->right;
        } else {
            section = branch->right;
            other = branch->left;
        }

        //

        nbh_s = neighborhood_(section, pt, rad, (level + 1) % dim);
        nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());
        if (dx2 < r2) {
            nbh_o = neighborhood_(other, pt, rad, (level + 1) % dim);
            nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
        }

        return nbh;
    };

   public:
    // default caller
    pointVec neighborhood(  //
        const point_t &pt,  //
        const double &rad) {
        size_t level = 0;
        return neighborhood_(root, pt, rad, level);
    }
};
