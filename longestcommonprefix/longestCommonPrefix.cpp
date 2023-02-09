// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include "longestCommonPrefix.h"
#include "../src/fmindexDBG.h"

struct LCPnode {

    FMPos pos;

    std::vector<length_t> boundaries;

    friend std::ostream& operator<<(std::ostream& o, const LCPnode& t);

    LCPnode() {
    }

    LCPnode(const FMPos& pos, const std::vector<length_t>& boundaries) {
        this->pos = pos;
        this->boundaries = boundaries;
    }
};

std::ifstream::pos_type filesize(std::string filename) {
    std::ifstream in(filename.c_str(),
                     std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

/*
 * file contains 'N' characters
 */
bool hasN(std::string filename) {

    std::ifstream i(filename);

    char c;

    while (i.get(c)) {

        if (c == 'N')
            return true;
    }

    return false;
}

uint64_t node_size(sa_node s) {
    return s.last - s.first_TERM;
}

uint64_t node_size(std::pair<sa_node, sa_node> p) {
    return node_size(p.first) + node_size(p.second);
}

uint64_t node_size(sa_node_n s) {
    return s.last - s.first_TERM;
}

uint64_t node_size(std::pair<sa_node_n, sa_node_n> p) {
    return node_size(p.first) + node_size(p.second);
}

void print_node(sa_node n) {

    std::cout << "[" << n.first_TERM << ", " << n.first_A << ", " << n.first_C
              << ", " << n.first_G << ", " << n.first_T << ", " << n.last << "]"
              << std::endl;
}

void print_node(sa_node_n n) {

    std::cout << "[" << n.first_TERM << ", " << n.first_A << ", " << n.first_C
              << ", " << n.first_G << ", " << n.first_N << ", " << n.first_T
              << ", " << n.last << "]" << std::endl;
}

sa_node merge_nodes(sa_node a, sa_node b) {

    assert(a.depth == b.depth);

    return {a.first_TERM + b.first_TERM,
            a.first_A + b.first_A,
            a.first_C + b.first_C,
            a.first_G + b.first_G,
            a.first_T + b.first_T,
            a.last + b.last,
            a.depth};
}

sa_node_n merge_nodes(sa_node_n a, sa_node_n b) {

    assert(a.depth == b.depth);

    return {a.first_TERM + b.first_TERM,
            a.first_A + b.first_A,
            a.first_C + b.first_C,
            a.first_G + b.first_G,
            a.first_N + b.first_N,
            a.first_T + b.first_T,
            a.last + b.last,
            a.depth};
}

inline uint64_t range_length(range_t r) {
    assert(r.second >= r.first);
    return r.second - r.first;
}

inline uint64_t leaf_size(sa_leaf L) {
    return range_length(L.rn);
}

inline uint64_t leaf_size(std::pair<sa_leaf, sa_leaf> P) {
    return leaf_size(P.first) + leaf_size(P.second);
}

void print_nodes(p_node p) {

    print_node(p.A);
    print_node(p.C);
    print_node(p.G);
    print_node(p.T);
}

range_t child_TERM(sa_node x) {
    return {x.first_TERM, x.first_A};
}
range_t child_A(sa_node x) {
    return {x.first_A, x.first_C};
}
range_t child_C(sa_node x) {
    return {x.first_C, x.first_G};
}
range_t child_G(sa_node x) {
    return {x.first_G, x.first_T};
}
range_t child_T(sa_node x) {
    return {x.first_T, x.last};
}

range_t child_TERM(sa_node_n x) {
    return {x.first_TERM, x.first_A};
}
range_t child_A(sa_node_n x) {
    return {x.first_A, x.first_C};
}
range_t child_C(sa_node_n x) {
    return {x.first_C, x.first_G};
}
range_t child_G(sa_node_n x) {
    return {x.first_G, x.first_N};
}
range_t child_N(sa_node_n x) {
    return {x.first_N, x.first_T};
}
range_t child_T(sa_node_n x) {
    return {x.first_T, x.last};
}

uint8_t number_of_children(sa_node N) {

    return uint8_t(N.last > N.first_T) + uint8_t(N.first_T > N.first_G) +
           uint8_t(N.first_G > N.first_C) + uint8_t(N.first_C > N.first_A) +
           uint8_t(N.first_A > N.first_TERM);
}

uint8_t number_of_children(sa_node_n N) {

    return uint8_t(N.last > N.first_T) + uint8_t(N.first_T > N.first_N) +
           uint8_t(N.first_N > N.first_G) + uint8_t(N.first_G > N.first_C) +
           uint8_t(N.first_C > N.first_A) + uint8_t(N.first_A > N.first_TERM);
}

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(sa_node N1, sa_node N2) {

    return uint8_t((N1.last > N1.first_T) or (N2.last > N2.first_T)) +
           uint8_t((N1.first_T > N1.first_G) or (N2.first_T > N2.first_G)) +
           uint8_t((N1.first_G > N1.first_C) or (N2.first_G > N2.first_C)) +
           uint8_t((N1.first_C > N1.first_A) or (N2.first_C > N2.first_A)) +
           uint8_t((N1.first_A > N1.first_TERM) or
                   (N2.first_A > N2.first_TERM));
}

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(sa_node_n N1, sa_node_n N2) {

    return uint8_t((N1.last > N1.first_T) or (N2.last > N2.first_T)) +
           uint8_t((N1.first_T > N1.first_N) or (N2.first_T > N2.first_N)) +
           uint8_t((N1.first_N > N1.first_G) or (N2.first_N > N2.first_G)) +
           uint8_t((N1.first_G > N1.first_C) or (N2.first_G > N2.first_C)) +
           uint8_t((N1.first_C > N1.first_A) or (N2.first_C > N2.first_A)) +
           uint8_t((N1.first_A > N1.first_TERM) or
                   (N2.first_A > N2.first_TERM));
}

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(std::pair<sa_node, sa_node> P) {

    return number_of_children(P.first, P.second);
}

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(std::pair<sa_node_n, sa_node_n> P) {

    return number_of_children(P.first, P.second);
}

template <class positionClass>
void computeLCPPrezza(FMIndexDBG<positionClass>* index, Bitvec2& LCP,
                      bool progress, uint& k) {
    // Initialize the LCP array with one more element than the length of the
    // original text
    LCP = (index->getTextLength() + 1);

    // cout << "\nNow navigating suffix tree leaves of size >= 2 to compute "
    //         "internal LCP values."
    //      << endl;

    uint64_t m = 0;      // portion of text covered by visited leaves
    uint64_t leaves = 0; // number of visited leaves
    uint64_t max_stack = 0;
    uint64_t lcp_values = 1; // number of filled LCP values

    {
        std::vector<FMPos> TMP_LEAVES;
        // TMP_LEAVES.reserve(5);

        std::vector<FMPos> S;
        S.reserve(1000);

        SARangePair rp = SARangePair(Range(0, index->getCounts()[1]),
                                     Range(0, index->getCounts()[1]));
        FMPos root = FMPos(rp, 0);
        S.push_back(root);

        int last_perc = -1;
        int perc = 0;

        while (!S.empty()) {
            max_stack = S.size() > max_stack ? S.size() : max_stack;

            FMPos L = S.back();
            S.pop_back();
            leaves++;

            assert(L.getRanges().getRangeSA().getEnd() >
                   L.getRanges().getRangeSA().getBegin());

            for (uint64_t i = L.getRanges().getRangeSA().getBegin() + 1;
                 i < L.getRanges().getRangeSA().getEnd(); ++i) {

                // assert(LCP[i] == nil);

                length_t h = L.getDepth();

                if (h < k) {
                    LCP[i] = 0;
                } else if (h == k) {
                    LCP[i] = 1;
                } else {
                    LCP[i] = 2;
                }

                lcp_values++;
                m++;
            }

            m++;

            assert(m <= index->getTextLength());

            next_leaves(index, L, TMP_LEAVES, 2);

            for (int i = TMP_LEAVES.size() - 1; i >= 0; --i) {
                S.push_back(TMP_LEAVES[i]);
            }

            if (progress) {
                perc = (100 * lcp_values) / index->getTextLength();

                if (perc > last_perc) {

                    std::cout << "LCP: " << perc << "%."
                              << "\r";
                    std::cout.flush();
                    last_perc = perc;
                }
            }
        }
    }

    // cout << "Visited leaves cover " << m << "/" << index->getTextLength()
    //      << " input characters." << endl;
    // cout << "Computed " << lcp_values << "/" << index->getTextLength()
    //      << " LCP values." << endl;

    // cout << "Max stack size = " << max_stack << endl;
    // cout << "Processed " << leaves << " suffix-tree leaves of size >= 2."
    //      << endl;

    // cout
    //     << "\nNow navigating suffix tree nodes to compute remaining LCP
    //     values."
    //     << endl;

    {
        std::vector<LCPnode> TMP_NODES;
        // TMP_NODES.reserve(5);

        uint64_t nodes = 0; // visited ST nodes
        max_stack = 0;

        std::vector<LCPnode> S;
        S.reserve(1000);

        LCPnode root;
        SARangePair rp = SARangePair(Range(0, index->getTextLength()),
                                     Range(0, index->getTextLength()));
        root.pos = FMPos(rp, 0);
        for (size_t i = 1; i < index->getSigma().size(); i++) {
            root.boundaries.emplace_back(index->getCounts()[i]);
        }
        S.push_back(root);

        int last_perc = -1;
        int perc = 0;

        while (!S.empty()) {

            max_stack = S.size() > max_stack ? S.size() : max_stack;

            LCPnode N = S.back();
            S.pop_back();
            nodes++;

            update_lcp(index, N, LCP, lcp_values, k);

            next_nodes(index, N, TMP_NODES);

            for (int i = TMP_NODES.size() - 1; i >= 0; --i) {
                S.push_back(TMP_NODES[i]);
            }

            TMP_NODES.clear();

            if (progress) {
                perc = (100 * lcp_values) / index->getTextLength();

                if (perc > last_perc) {

                    std::cout << "LCP: " << perc << "%."
                              << "\r";
                    std::cout.flush();
                    last_perc = perc;
                }
            }
        }
    }
    if (progress) {

        std::cout << "LCP: 100%."
                  << "\n";
        std::cout << "Maximum stack size: " << max_stack << std::endl;
    }
}

template <class positionClass>
void next_leaves(FMIndexDBG<positionClass>* index, FMPos& L,
                 std::vector<FMPos>& TMP_LEAVES, length_t min_n_children) {

    for (length_t i = 1; i < index->getSigma().size(); i++) {
        SARangePair newRanges;
        index->findRangesWithExtraCharBackward(i, L.getRanges(), newRanges);
        if (newRanges.getRangeSA().width() >= min_n_children) {
            TMP_LEAVES.emplace_back(newRanges, L.getDepth() + 1);
        }
    }

    std::sort(TMP_LEAVES.begin(), TMP_LEAVES.end(),
              [](const FMPos& lhs, const FMPos& rhs) {
                  return lhs.getRanges().getRangeSA().width() <
                         rhs.getRanges().getRangeSA().width();
              });
}

template <class positionClass>
void update_lcp(FMIndexDBG<positionClass>* index, LCPnode& x, Bitvec2& LCP,
                uint64_t& lcp_values, uint& k) {
    for (size_t i = 0; i < index->getSigma().size() - 1; i++) {
        if ((x.boundaries[i] > ((i == 0)
                                    ? x.pos.getRanges().getRangeSA().getBegin()
                                    : x.boundaries[i - 1])) &&
            (x.boundaries[i] != x.pos.getRanges().getRangeSA().getEnd())) {
            length_t h = x.pos.getDepth();

            if (h < k) {
                LCP[x.boundaries[i]] = 0;
            } else if (h == k) {
                LCP[x.boundaries[i]] = 1;
            } else {
                LCP[x.boundaries[i]] = 2;
            }

            lcp_values++;
        }
    }
}

template <class positionClass>
void next_nodes(FMIndexDBG<positionClass>* index, LCPnode& x,
                std::vector<LCPnode>& TMP_NODES) {
    for (size_t i = 1; i < index->getSigma().size(); i++) {
        SARangePair newRanges;
        index->findRangesWithExtraCharBackward(i, x.pos.getRanges(), newRanges);
        FMPos newPos(newRanges, x.pos.getDepth() + 1);
        // LCPnode result;
        // result.pos = newPos;
        size_t childrenCounter = 0;
        std::vector<length_t> boundaries(5, 0);
        for (size_t j = 1; j < index->getSigma().size(); j++) {
            childrenCounter += index->findRangesWithExtraCharForward(
                j, newPos.getRanges(), newRanges);
            boundaries[j - 1] = newRanges.getRangeSA().getBegin();
        }
        childrenCounter +=
            newPos.getRanges().getRangeSA().getBegin() != boundaries[0];
        if (childrenCounter >= 2) {
            TMP_NODES.emplace_back(newPos, boundaries);
        }
    }

    std::sort(TMP_NODES.begin(), TMP_NODES.end(),
              [](const LCPnode& lhs, const LCPnode& rhs) {
                  return lhs.pos.getRanges().width() <
                         rhs.pos.getRanges().width();
              });

    // for (size_t i = 1; i < index->getSigma().size(); i++) {
    //     Range r = Range(x.front().getRanges().getRangeSA().getBegin(),
    //                     x.back().getRanges().getRangeSA().getEnd());
    //     SARangePair rp = SARangePair(r, r);
    //     FMPos currentPos = x[i];
    //     size_t childrenCounter = 0;
    //     toBeNamed childVector;
    //     for (size_t j = 0; j < index->getSigma().size(); j++) {
    //         SARangePair newRanges;
    //         childrenCounter += index->findRangesWithExtraCharForward(
    //             j, currentPos.getRanges(), newRanges);
    //         childVector.emplace_back(newRanges, currentPos.getDepth() + 1);
    //     }
    //     if (childrenCounter >= 2) {
    //         TMP_NODES.push_back(childVector);
    //     }
    // }

    // std::sort(TMP_NODES.begin(), TMP_NODES.end(),
    //           [](const toBeNamed& lhs, const toBeNamed& rhs) {
    //               return lhs.back().getRanges().getRangeSA().getEnd() -
    //                          lhs.front().getRanges().getRangeSA().getBegin()
    //                          <
    //                      rhs.back().getRanges().getRangeSA().getEnd() -
    //                          rhs.front().getRanges().getRangeSA().getBegin();
    //           });

    // auto stub = 0;
}

template void computeLCPPrezza(FMIndexDBG<FMPos>* index, Bitvec2& LCP,
                               bool progress, uint& k);
template void computeLCPPrezza(FMIndexDBG<FMPosSFR>* index, Bitvec2& LCP,
                               bool progress, uint& k);