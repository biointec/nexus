// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#ifndef INCLUDE_HPP_
#define INCLUDE_HPP_

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

struct LCPnode;
template <class positionClass> class FMIndexDBG;
class Bitvec2;
class FMPos;

typedef std::pair<uint64_t, uint64_t> range_t;
typedef uint64_t length_t;

std::ifstream::pos_type filesize(std::string filename);

/*
 * representation of a right-maximal substring (SA node) as a list of BWT
 * intervals
 */
struct sa_node {

    // right-maximal substring: string W such that Wa_1, ..., Wa_k occur in the
    // text for at least k>=2 characters a_1, ..., a_k

    uint64_t first_TERM;
    uint64_t first_A;
    uint64_t first_C;
    uint64_t first_G;
    uint64_t first_T;
    uint64_t last;

    // depth = |W|
    uint64_t depth;

    uint64_t key() {
        return first_TERM;
    }
};

struct sa_node_n {

    // right-maximal substring: string W such that Wa_1, ..., Wa_k occur in the
    // text for at least k>=2 characters a_1, ..., a_k

    uint64_t first_TERM;
    uint64_t first_A;
    uint64_t first_C;
    uint64_t first_G;
    uint64_t first_N;
    uint64_t first_T;
    uint64_t last;

    // depth = |W|
    uint64_t depth;

    uint64_t key() {
        return first_TERM;
    }
};

/*
 * file contains 'N' characters
 */
bool hasN(std::string filename);

uint64_t node_size(sa_node s);

uint64_t node_size(std::pair<sa_node, sa_node> p);

uint64_t node_size(sa_node_n s);

uint64_t node_size(std::pair<sa_node_n, sa_node_n> p);

void print_node(sa_node n);

void print_node(sa_node_n n);

sa_node merge_nodes(sa_node a, sa_node b);

sa_node_n merge_nodes(sa_node_n a, sa_node_n b);

/*
 * suffix array leaf = BWT range (inclusive) of W.TERM, for some string W.
 *
 */
struct sa_leaf {

    // rn.first = first position of range. Equivalently, number of suffixes
    // smaller than W.TERM (valid also if W.TERM does not occur) rn.second =
    // last position (excluded) of interval.  Equivalently, number of suffixes
    // smaller than W.TERM + number of occurrences of W.TERM if last == first,
    // then W.TERM does not occur (however, 'first' is in any case number of
    // suffixes smaller than W.TERM)
    range_t rn;

    // depth = |W.TERM|
    uint64_t depth;

    uint64_t key() {
        return rn.first;
    }
};

inline uint64_t range_length(range_t r);

inline uint64_t leaf_size(sa_leaf L);

inline uint64_t leaf_size(std::pair<sa_leaf, sa_leaf> P);

struct p_range {

    range_t A;
    range_t C;
    range_t G;
    range_t T;
};

struct p_node {

    sa_node A;
    sa_node C;
    sa_node G;
    sa_node T;
};

struct p_range_n {

    range_t A;
    range_t C;
    range_t G;
    range_t N;
    range_t T;
};

struct p_node_n {

    sa_node_n A;
    sa_node_n C;
    sa_node_n G;
    sa_node_n N;
    sa_node_n T;
};

void print_nodes(p_node p);

struct p_rank {

  public:
    uint64_t A;
    uint64_t C;
    uint64_t G;
    uint64_t T;

    p_rank operator+(const p_rank& a) const {

        return {a.A + A, a.C + C, a.G + G, a.T + T};
    }

    bool operator==(const p_rank& a) const {

        return a.A == A and a.C == C and a.G == G and a.T == T;
    }

    bool operator!=(const p_rank& a) const {

        return a.A != A or a.C != C or a.G != G or a.T != T;
    }

    bool operator<=(const p_rank& a) const {

        return A <= a.A and C <= a.C and G <= a.G and T <= a.T;
    }
};

struct p_rank_n {

  public:
    uint64_t A;
    uint64_t C;
    uint64_t G;
    uint64_t N;
    uint64_t T;

    p_rank_n operator+(const p_rank_n& a) const {

        return {a.A + A, a.C + C, a.G + G, a.N + N, a.T + T};
    }

    bool operator==(const p_rank_n& a) const {

        return a.A == A and a.C == C and a.G == G and a.N == N and a.T == T;
    }

    bool operator!=(const p_rank_n& a) const {

        return a.A != A or a.C != C or a.G != G or a.N != N or a.T != T;
    }

    bool operator<=(const p_rank_n& a) const {

        return A <= a.A and C <= a.C and G <= a.G and N <= a.N and T <= a.T;
    }
};

inline p_range fold_ranks(p_rank& a, p_rank& b) {

    return {{a.A, b.A}, {a.C, b.C}, {a.G, b.G}, {a.T, b.T}};
}

inline p_range_n fold_ranks(p_rank_n& a, p_rank_n& b) {

    return {{a.A, b.A}, {a.C, b.C}, {a.G, b.G}, {a.N, b.N}, {a.T, b.T}};
}

inline uint64_t popcount128(__uint128_t x) {

    return __builtin_popcountll(uint64_t(x >> 64)) +
           __builtin_popcountll(x & 0xFFFFFFFFFFFFFFFF);
}

range_t child_TERM(sa_node x);
range_t child_A(sa_node x);
range_t child_C(sa_node x);
range_t child_G(sa_node x);
range_t child_T(sa_node x);

range_t child_TERM(sa_node_n x);
range_t child_A(sa_node_n x);
range_t child_C(sa_node_n x);
range_t child_G(sa_node_n x);
range_t child_N(sa_node_n x);
range_t child_T(sa_node_n x);

inline bool has_child_TERM(sa_node N) {
    return N.first_A > N.first_TERM;
}
inline bool has_child_A(sa_node N) {
    return N.first_C > N.first_A;
}
inline bool has_child_C(sa_node N) {
    return N.first_G > N.first_C;
}
inline bool has_child_G(sa_node N) {
    return N.first_T > N.first_G;
}
inline bool has_child_T(sa_node N) {
    return N.last > N.first_T;
}

inline bool has_child_TERM(sa_node_n N) {
    return N.first_A > N.first_TERM;
}
inline bool has_child_A(sa_node_n N) {
    return N.first_C > N.first_A;
}
inline bool has_child_C(sa_node_n N) {
    return N.first_G > N.first_C;
}
inline bool has_child_G(sa_node_n N) {
    return N.first_N > N.first_G;
}
inline bool has_child_N(sa_node_n N) {
    return N.first_T > N.first_N;
}
inline bool has_child_T(sa_node_n N) {
    return N.last > N.first_T;
}

uint8_t number_of_children(sa_node N);

uint8_t number_of_children(sa_node_n N);

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(sa_node N1, sa_node N2);

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(sa_node_n N1, sa_node_n N2);

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(std::pair<sa_node, sa_node> P);

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(std::pair<sa_node_n, sa_node_n> P);

template <typename lcp_int_t>
void update_lcp(sa_node& x, std::vector<lcp_int_t>& LCP, uint64_t& lcp_values) {

    assert(x.first_A >= x.first_TERM);
    assert(x.first_C >= x.first_A);
    assert(x.first_G >= x.first_C);
    assert(x.first_T >= x.first_G);

    assert(number_of_children(x) >= 2);

    lcp_int_t nil = ~lcp_int_t(0);

    if (has_child_TERM(x) and x.first_A != x.last) {
        assert(LCP[x.first_A] == nil);
        LCP[x.first_A] = x.depth;
        lcp_values++;
    }
    if (has_child_A(x) and x.first_C != x.last) {
        assert(LCP[x.first_C] == nil);
        LCP[x.first_C] = x.depth;
        lcp_values++;
    }
    if (has_child_C(x) and x.first_G != x.last) {
        assert(LCP[x.first_G] == nil);
        LCP[x.first_G] = x.depth;
        lcp_values++;
    }
    if (has_child_G(x) and x.first_T != x.last) {
        assert(LCP[x.first_T] == nil);
        LCP[x.first_T] = x.depth;
        lcp_values++;
    }
}

template <typename lcp_int_t>
void update_lcp(sa_node_n& x, std::vector<lcp_int_t>& LCP,
                uint64_t& lcp_values) {

    assert(x.first_A >= x.first_TERM);
    assert(x.first_C >= x.first_A);
    assert(x.first_G >= x.first_C);
    assert(x.first_N >= x.first_G);
    assert(x.first_T >= x.first_N);

    assert(number_of_children(x) >= 2);

    lcp_int_t nil = ~lcp_int_t(0);

    if (has_child_TERM(x) and x.first_A != x.last) {
        assert(LCP[x.first_A] == nil);
        LCP[x.first_A] = x.depth;
        lcp_values++;
    }
    if (has_child_A(x) and x.first_C != x.last) {
        assert(LCP[x.first_C] == nil);
        LCP[x.first_C] = x.depth;
        lcp_values++;
    }
    if (has_child_C(x) and x.first_G != x.last) {
        assert(LCP[x.first_G] == nil);
        LCP[x.first_G] = x.depth;
        lcp_values++;
    }
    if (has_child_G(x) and x.first_N != x.last) {
        assert(LCP[x.first_N] == nil);
        LCP[x.first_N] = x.depth;
        lcp_values++;
    }
    if (has_child_N(x) and x.first_T != x.last) {
        assert(LCP[x.first_T] == nil);
        LCP[x.first_T] = x.depth;
        lcp_values++;
    }
}

template <class positionClass>
void computeLCPPrezza(FMIndexDBG<positionClass>* index, Bitvec2& LCP,
                      bool progress, uint& k);

template <class positionClass>
void next_leaves(FMIndexDBG<positionClass>* index, FMPos& L,
                 std::vector<FMPos>& TMP_LEAVES, length_t min_n_children);

template <class positionClass>
void update_lcp(FMIndexDBG<positionClass>* index, LCPnode& x, Bitvec2& LCP,
                uint64_t& lcp_values, uint& k);

template <class positionClass>
void next_nodes(FMIndexDBG<positionClass>* index, LCPnode& x,
                std::vector<LCPnode>& TMP_NODES);

#endif /* INCLUDE_HPP_ */
