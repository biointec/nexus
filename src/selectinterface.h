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

#ifndef SELECT_INTERFACE_H
#define SELECT_INTERFACE_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

#include "../sux/bits/EliasFano.hpp"
#include "../sux/bits/Rank9Sel.hpp"
#include "../sux/bits/Select.hpp"
#include "../sux/bits/SimpleSelect.hpp"
#include "../sux/bits/SimpleSelectHalf.hpp"
#pragma GCC diagnostic pop

#include "bitvec.h"

enum class SelectOption { ELIASFANO, SELECT9, SIMPLE, HALF };

inline const std::string tostring(SelectOption v) {
    switch (v) {
    case SelectOption::ELIASFANO:
        return "ELIASFANO";
    case SelectOption::SELECT9:
        return "SELECT9";
    case SelectOption::SIMPLE:
        return "SIMPLE";
    case SelectOption::HALF:
        return "HALF";

    default:
        return "NOT AN OPTION!!";
    }
}

class SelectInterface {
  private:
    sux::Select* selectClass = nullptr;
    SelectOption option;
    int max_log2_longwords_per_subinventory = 0;

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
        delete selectClass;

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
     * Swap two SelectInterface objects
     * @param lhs Left hand size
     * @param rhs Right hand size
     */
    friend void swap(SelectInterface& lhs, SelectInterface& rhs) {
        using std::swap;
        swap(lhs.N, rhs.N);
        swap(lhs.numWords, rhs.numWords);
        swap(lhs.bv, rhs.bv);
    }

  public:
    SelectInterface(SelectOption option = SelectOption::SIMPLE, uint64_t N = 0,
                    const int max_log2_longwords_per_subinventory = 0)
        : option(option), max_log2_longwords_per_subinventory(
                              max_log2_longwords_per_subinventory),
          N(N), bv(NULL) {
        allocateMem();
    }

    /**
     * Move constructor
     * @param rhs Right hand size
     */
    SelectInterface(SelectInterface&& rhs) : SelectInterface() {
        swap(*this, rhs);
    }

    /**
     * Move assignment operator (shallow copy)
     * @param rhs Right hand size
     */
    SelectInterface& operator=(SelectInterface&& rhs) {
        swap(*this, rhs);
        return *this;
    }

    /**
     * Deleted copy constructor and copy assignment operator
     */
    SelectInterface(const SelectInterface&) = delete;
    SelectInterface& operator=(const SelectInterface&) = delete;

    ~SelectInterface() {
        delete selectClass;
        free(bv);
    }

    void setN(uint64_t newN) {
        N = newN;
        allocateMem();
    }

    void setOption(const SelectOption& newOption) {
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
        // TODO other options
        switch (option) {
        case SelectOption::SIMPLE:
            selectClass = new sux::bits::SimpleSelect<>(
                bv, N, max_log2_longwords_per_subinventory);
            break;
        case SelectOption::HALF:
            selectClass = new sux::bits::SimpleSelectHalf<>(bv, N);
            break;
        case SelectOption::SELECT9:
            selectClass = new sux::bits::Rank9Sel<>(bv, N);
            break;
        case SelectOption::ELIASFANO:
            selectClass = new sux::bits::EliasFano<>(bv, N);
            break;
        default:
            throw std::runtime_error("SelectOption " + tostring(option) +
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

    std::size_t select(uint64_t rank) const {
        return selectClass->select(rank);
    }

    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv, numWords * sizeof(uint64_t));
        selectClass->write(ofs);
    }

    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));
        allocateMem();
        ifs.read((char*)bv, numWords * sizeof(uint64_t));
        switch (option) {
        case SelectOption::SIMPLE:
            selectClass = new sux::bits::SimpleSelect<>(ifs, bv);
            break;
        case SelectOption::HALF:
            selectClass = new sux::bits::SimpleSelectHalf<>(ifs, bv);
            break;
        case SelectOption::SELECT9:
            selectClass = new sux::bits::Rank9Sel<>(ifs, bv, N);
            break;
        case SelectOption::ELIASFANO:
            selectClass = new sux::bits::EliasFano<>(ifs, N);
            break;
        default:
            throw std::runtime_error("SelectOption " + tostring(option) +
                                     " is not supported yet!");
        }
    }
};

#endif