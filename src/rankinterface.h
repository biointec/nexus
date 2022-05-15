/******************************************************************************
 *  Nexus: Pan-genome compacted de Bruijn graphs with support for approximate *
 *      pattern matching using search schemes                                 *
 *                                                                            *
 *  Copyright (C) 2022 - Lore Depuydt <lore.depuydt@ugent.be>,                *
 *                       Luca Renders <luca.renders@ugent.be> and             *
 *                       Jan Fostier <jan.fostier@ugent.be>                   *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#ifndef RANK_INTERFACE_H
#define RANK_INTERFACE_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include "../sux/bits/EliasFano.hpp"
#include "../sux/bits/Rank.hpp"
#include "../sux/bits/Rank9.hpp"
#include "../sux/bits/Rank9Sel.hpp"
#pragma GCC diagnostic pop

#include "bitvec.h"

enum class RankOption { ELIASFANO, RANK9, RANK9SELECT };

inline const std::string tostring(RankOption v) {
    switch (v) {
    case RankOption::ELIASFANO:
        return "ELIASFANO";
    case RankOption::RANK9:
        return "RANK9";
    case RankOption::RANK9SELECT:
        return "RANK9SELECT";

    default:
        return "NOT AN OPTION!!";
    }
}

class RankInterface {
  private:
    sux::Rank* rankClass = nullptr;
    RankOption option;

    uint64_t N;        // size of the bitvector
    uint64_t numWords; // number of words in the bitvector
    uint64_t* bv;      // interleaved bitvectors

    /**
     * Allocate memory for bv
     */
    void allocateMem() {
        // free existing allocations
        free(bv);
        bv = NULL;
        delete rankClass;

        if (N == 0) { // special case for N == 0
            numWords = 0;
            return;
        }

        const size_t B = 64; // memory alignment in bytes

        // allocate memory for the bitvector (and set to zero)
        numWords = ((N + 63) / 64);
        // numBytes must be an integral multiple of B
        uint64_t numBytes = ((numWords * sizeof(uint64_t) + B - 1) / B) * B;
        bv = (uint64_t*)aligned_alloc(B, numBytes);
        memset((void*)bv, 0, numBytes);
    }

    /**
     * Swap two RankInterface objects
     * @param lhs Left hand size
     * @param rhs Right hand size
     */
    friend void swap(RankInterface& lhs, RankInterface& rhs) {
        using std::swap;
        swap(lhs.N, rhs.N);
        swap(lhs.numWords, rhs.numWords);
        swap(lhs.bv, rhs.bv);
    }

  public:
    RankInterface(RankOption option = RankOption::RANK9, uint64_t N = 0)
        : option(option), N(N), bv(NULL) {
        allocateMem();
    }

    /**
     * Move constructor
     * @param rhs Right hand size
     */
    RankInterface(RankInterface&& rhs) : RankInterface() {
        swap(*this, rhs);
    }

    /**
     * Move assignment operator (shallow copy)
     * @param rhs Right hand size
     */
    RankInterface& operator=(RankInterface&& rhs) {
        swap(*this, rhs);
        return *this;
    }

    /**
     * Deleted copy constructor and copy assignment operator
     */
    RankInterface(const RankInterface&) = delete;
    RankInterface& operator=(const RankInterface&) = delete;

    ~RankInterface() {
        delete rankClass;
        free(bv);
    }

    void setN(uint64_t newN) {
        N = newN;
        allocateMem();
    }

    void setOption(const RankOption& newOption) {
        option = newOption;
    }

    /**
     * Get a bit at position p
     * @param p Position
     * @return true or false
     */
    bool operator[](size_t p) const {

        assert(p < N);
        size_t w = (p / 64);
        size_t b = p % 64;
        return (bv[w] & (1ull << b)) != 0;
    }

    /**
     * Get a bit reference at position p
     * @param p Position
     * @return Bit reference object
     */
    Bitref operator[](size_t p) {

        assert(p < N);
        size_t w = (p / 64);
        size_t b = p % 64;
        return Bitref(bv[w], 1ull << b);
    }

    void indexInterface() {
        // create the appropriate Select datastructure
        switch (option) {
        case RankOption::RANK9:
            rankClass = new sux::bits::Rank9<>(bv, N);
            break;
        case RankOption::RANK9SELECT:
            rankClass = new sux::bits::Rank9Sel<>(bv, N);
            break;
        case RankOption::ELIASFANO:
            rankClass = new sux::bits::EliasFano<>(bv, N);
            break;
        default:
            throw std::runtime_error("RankOption " + tostring(option) +
                                     " is not supported yet!");
        }
    }

    /**
     * Return the size of the underlying bitvector
     * @return The size of the bitvector
     */
    size_t size() const {
        return N;
    }

    uint64_t rank(const std::size_t k) const {
        return rankClass->rank(k);
    }

    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv, numWords * sizeof(uint64_t));
        rankClass->write(ofs);
    }

    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));
        allocateMem();
        ifs.read((char*)bv, numWords * sizeof(uint64_t));
        switch (option) {
        case RankOption::RANK9:
            rankClass = new sux::bits::Rank9<>(ifs, bv, N);
            break;
        case RankOption::RANK9SELECT:
            rankClass = new sux::bits::Rank9Sel<>(ifs, bv, N);
            break;
        case RankOption::ELIASFANO:
            rankClass = new sux::bits::EliasFano<>(ifs, N);
            break;
        default:
            throw std::runtime_error("RankOption " + tostring(option) +
                                     " is not supported yet!");
        }
    }
};

#endif