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

#ifndef BITVEC_H
#define BITVEC_H

/**
 * The implementation implements the rank9 algorithm as described in
 * S. Vigna, "Broadword Implementation of Rank/Select Queries", WEA 2008
 * It relies on GCC's __builtin_popcountll, so please build this software
 * using the -mpopcnt flag to enable the SSE 4.2 POPCNT instruction.
 */

#include "alphabet.h"

#include <cassert>
#include <cmath> // taking the log
#include <cstdint>
#include <cstring>
#include <fstream> // I/O
#include <vector>

// ============================================================================
// BIT REFERENCE CLASS
// ============================================================================

class Bitref {

  private:
    size_t& wordRef; // reference to a word in the bitvector
    size_t bitmask;  // bitmask of the form (1 << bitIdx)

  public:
    /**
     * Constructor
     * @param wordRef Reference to a word in the bitvector
     * @param bitmask Bitmask of the form (1 << bitIdx)
     */
    Bitref(size_t& wordRef, size_t bitmask)
        : wordRef(wordRef), bitmask(bitmask) {
    }

    /**
     * Set bitref to particular value
     * @param val Target value
     * @return Bitref reference after modification
     */
    const Bitref& operator=(bool val) {
        if (val)
            wordRef |= bitmask;
        else
            wordRef &= ~bitmask;
        return *this;
    }

    /**
     * Set bitref to another bitref
     * @param br Another bitref
     * @return Bitref reference
     */
    const Bitref& operator=(const Bitref& br) {
        return this->operator=(bool(br));
    }

    /**
     * Bool conversion operator
     */
    operator bool() const {
        return (wordRef & bitmask) != 0;
    }
};

// ============================================================================
// BIT VECTOR CLASS
// ============================================================================

class Bitvec {

  private:
    size_t N;                   // size of the bitvector
    std::vector<size_t> bv;     // actual bitvector
    std::vector<size_t> counts; // interleaved 1st and 2nd level counts

  public:
    /**
     * Get a bit at a certain position
     * @param p Position
     * @return true or false
     */
    bool operator[](size_t p) const {
        assert(p < N);
        size_t w = p / 64;
        size_t b = p % 64;
        return (bv[w] & (1ull << b)) != 0;
    }

    /**
     * Get a bit reference at a certain position
     * @param p Position
     * @return Bit reference object
     */
    Bitref operator[](size_t p) {
        assert(p < N);
        size_t w = p / 64;
        size_t b = p % 64;
        return Bitref(bv[w], 1ull << b);
    }

    /**
     * Create an index for the bitvector to support fast rank operations
     */
    void index() {
        counts = std::vector<size_t>((bv.size() + 7) / 4, 0ull);

        size_t countL1 = 0, countL2 = 0;
        for (size_t w = 0, q = 0; w < bv.size(); w++) {
            if (w % 8 == 0) { // store the L1 counts
                countL1 += countL2;
                counts[q] = countL1;
                countL2 = __builtin_popcountll(bv[w]);
                q += 2;
            } else { // store the L2 counts
                counts[q - 1] |= (countL2 << (((w % 8) - 1) * 9));
                countL2 += __builtin_popcountll(bv[w]);
            }
        }
    }

    /**
     * Get the number of 1-bits within the range [0...p[ (preceding pos p)
     * @param p Position
     */
    size_t rank(size_t p) const {
        assert(p < N);
        size_t w = p / 64;      // word index
        size_t b = p % 64;      // bit offset
        size_t q = (w / 8) * 2; // counts index

        // add the first-level counts
        size_t rv = counts[q];

        // add the second-level counts
        int64_t t = (w % 8) - 1;
        rv += counts[q + 1] >> (t + (t >> 60 & 8)) * 9 & 0x1FF;

        // add the popcount in the final word
        return rv + __builtin_popcountll((bv[w] << 1) << (63 - b));
    }

    /**
     * Write the bitvector to an open filestream
     * @param ofs Open output filestream
     */
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv.data(), bv.size() * sizeof(size_t));
        ofs.write((char*)counts.data(), counts.size() * sizeof(size_t));
    }

    /**
     * Read the bitvector from an open filestream
     * @param ifs Open input filestream
     */
    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));

        bv.resize((N + 63) / 64);
        ifs.read((char*)bv.data(), bv.size() * sizeof(size_t));

        counts.resize((bv.size() + 7) / 4);
        ifs.read((char*)counts.data(), counts.size() * sizeof(size_t));
    }

    /**
     * Return the size of the bitvector
     * @return The size of the bitvector
     */
    size_t size() const {
        return N;
    }

    /**
     * Default constructor, move constructor and move assignment operator
     */
    Bitvec() : N(0){};
    Bitvec(Bitvec&& rhs) = default;
    Bitvec& operator=(Bitvec&& rhs) = default;

    /**
     * Deleted copy constructor and copy assignment operator
     */
    Bitvec(const Bitvec&) = delete;
    Bitvec& operator=(const Bitvec&) = delete;

    /**
     * Constructor
     * @param N Number of bits in the bitvector
     */
    Bitvec(size_t N) : N(N), bv((N + 63) / 64, 0ull) {
    }
};

// ============================================================================
// TEMPLATED INTERLEAVED BIT VECTOR CLASS
// ============================================================================

template <size_t S> // S is the size of the alphabet
class BitvecIntl {

  private:
    size_t N;          // size of the bitvector
    size_t bvSize;     // number of words in the bitvector
    size_t* bv;        // interleaved bitvectors
    size_t countsSize; // number of words in the counts vector
    size_t* counts;    // interleaved 1st and 2nd level counts

    /**
     * Allocate memory for bv and counts
     */
    void allocateMem() {
        // free existing allocations
        free(bv);
        bv = NULL;
        free(counts);
        counts = NULL;

        if (N == 0) { // special case for N == 0
            bvSize = countsSize = 0;
            return;
        }

        const size_t B = 64; // memory alignment in bytes

        // allocate memory for the bitvector (and set to zero)
        bvSize = S * ((N + 63) / 64);
        // numBytes must be an integral multiple of B
        size_t numBytes = ((bvSize * sizeof(size_t) + B - 1) / B) * B;
        bv = (size_t*)aligned_alloc(B, numBytes);
        memset((void*)bv, 0, numBytes);

        // allocate memory for the counts
        countsSize = 2 * S * ((N + 511) / 512);
        // numBytes must be an integral multiple of B
        numBytes = ((countsSize * sizeof(size_t) + B - 1) / B) * B;
        counts = (size_t*)aligned_alloc(B, numBytes);
    }

    /**
     * Swap two BitvecIntl objects
     * @param lhs Left hand size
     * @param rhs Right hand size
     */
    friend void swap(BitvecIntl<S>& lhs, BitvecIntl<S>& rhs) {
        using std::swap;
        swap(lhs.N, rhs.N);
        swap(lhs.bvSize, rhs.bvSize);
        swap(lhs.bv, rhs.bv);
        swap(lhs.countsSize, rhs.countsSize);
        swap(lhs.counts, rhs.counts);
    }

  public:
    /**
     * Get a bit at position p for character c
     * @param c Character index [0,1,...,S[
     * @param p Position
     * @return true or false
     */
    bool operator()(size_t c, size_t p) const {
        assert(c < S);
        assert(p < N);
        size_t w = (p / 64) * S + c;
        size_t b = p % 64;
        return (bv[w] & (1ull << b)) != 0;
    }

    /**
     * Get a bit reference at position p for character c
     * @param c Character index [0,1,...,S[
     * @param p Position
     * @return Bit reference object
     */
    Bitref operator()(size_t c, size_t p) {
        assert(c < S);
        assert(p < N);
        size_t w = (p / 64) * S + c;
        size_t b = p % 64;
        return Bitref(bv[w], 1ull << b);
    }

    /**
     * Create an index for the bitvector to support fast rank operations
     */
    void index() {
        // reset counts to zero
        memset((void*)counts, 0, countsSize * sizeof(size_t));

        for (size_t c = 0; c < S; c++) {
            size_t countL1 = 0, countL2 = 0;
            for (size_t w = c, q = 2 * c; w < bvSize; w += S) {
                size_t numBits = __builtin_popcountll(bv[w]);
                if (w % (8 * S) == c) { // store the L1 counts
                    countL1 += countL2;
                    counts[q] = countL1;
                    countL2 = numBits;
                    q += 2 * S;
                } else { // store the L2 counts
                    size_t L2offs = 9 * ((w / S % 8) - 1);
                    counts[q + 1 - 2 * S] |= (countL2 << L2offs);
                    countL2 += numBits;
                }
            }
        }
    }

    /**
     * Get the number of 1-bits within the range [0...p[ (preceding pos p)
     * @param c Character index [0,1,...,S[
     * @param p Position
     */
    size_t rank(size_t c, size_t p) const {
        assert(c < S);
        assert(p < N);
        size_t w = (p / 64) * S + c;          // word index
        size_t b = p % 64;                    // bit offset
        size_t q = (p / 512) * 2 * S + 2 * c; // counts index

        // add the first-level counts
        size_t rv = counts[q];

        // add the second-level counts
        int64_t t = ((p / 64) % 8) - 1;
        rv += counts[q + 1] >> (t + (t >> 60 & 8)) * 9 & 0x1FF;

        // add the popcount in the final word
        return rv + __builtin_popcountll((bv[w] << 1) << (63 - b));
    }

    /**
     * Write the bitvector to an open filestream
     * @param ofs Open output filestream
     */
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv, bvSize * sizeof(size_t));
        ofs.write((char*)counts, countsSize * sizeof(size_t));
    }

    /**
     * Read the bitvector from an open filestream
     * @param ifs Open input filestream
     */
    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));
        allocateMem();
        ifs.read((char*)bv, bvSize * sizeof(size_t));
        ifs.read((char*)counts, countsSize * sizeof(size_t));
    }

    /**
     * Return the size of the bitvector
     * @return The size of the bitvector
     */
    size_t size() const {
        return N;
    }

    /**
     * Constructor
     * @param N Number of bits in the interleaved bitvector per character
     */
    BitvecIntl(size_t N = 0) : N(N), bv(NULL), counts(NULL) {
        allocateMem();
    }

    /**
     * Move constructor
     * @param rhs Right hand size
     */
    BitvecIntl(BitvecIntl<S>&& rhs) : BitvecIntl() {
        swap(*this, rhs);
    }

    /**
     * Move assignment operator (shallow copy)
     * @param rhs Right hand size
     */
    BitvecIntl& operator=(BitvecIntl<S>&& rhs) {
        swap(*this, rhs);
        return *this;
    }

    /**
     * Deleted copy constructor and copy assignment operator
     */
    BitvecIntl(const BitvecIntl<S>&) = delete;
    BitvecIntl& operator=(const BitvecIntl<S>&) = delete;

    /**
     * Destructor
     */
    ~BitvecIntl() {
        free(bv);
        free(counts);
    }
};

// ============================================================================
// TEMPLATED NON-INTERLEAVED BIT VECTOR CLASS
// ============================================================================

template <size_t S> // S is the size of the alphabet
class BitvecNonIntl {

  protected:
    size_t N;          // size of the bitvector
    size_t bvSize;     // number of words in the bitvector
    size_t* bv;        // interleaved bitvectors
    size_t countsSize; // number of words in the counts vector
    size_t* counts;    // interleaved 1st and 2nd level counts

    /**
     * Allocate memory for bv and counts
     */
    void allocateMem() {
        // free existing allocations
        free(bv);
        bv = NULL;
        free(counts);
        counts = NULL;

        if (N == 0) { // special case for N == 0
            bvSize = countsSize = 0;
            return;
        }

        const size_t B = 64; // memory alignment in bytes

        // allocate memory for the bitvector (and set to zero)
        bvSize = S * ((N + 63) / 64);
        // numBytes must be an integral multiple of B
        size_t numBytes = ((bvSize * sizeof(size_t) + B - 1) / B) * B;
        bv = (size_t*)aligned_alloc(B, numBytes);
        memset((void*)bv, 0, numBytes);

        // allocate memory for the counts
        countsSize = S * (((N + 64 * 8 - 1) / (64 * 8)) * 2 + 2);
        // numBytes must be an integral multiple of B
        numBytes = ((countsSize * sizeof(size_t) + B - 1) / B) * B;
        counts = (size_t*)aligned_alloc(B, numBytes);
    }

    /**
     * Swap two BitvecIntl objects
     * @param lhs Left hand size
     * @param rhs Right hand size
     */
    friend void swap(BitvecNonIntl<S>& lhs, BitvecNonIntl<S>& rhs) {
        using std::swap;
        swap(lhs.N, rhs.N);
        swap(lhs.bvSize, rhs.bvSize);
        swap(lhs.bv, rhs.bv);
        swap(lhs.countsSize, rhs.countsSize);
        swap(lhs.counts, rhs.counts);
    }

  public:
    /**
     * Get a bit at position p for character c
     * @param c Character index [0,1,...,S[
     * @param p Position
     * @return true or false
     */
    bool operator()(size_t c, size_t p) const {
        assert(c < S);
        assert(p < N);
        size_t w = (bvSize / S) * c + (p / 64);
        size_t b = p % 64;
        return (bv[w] & (1ull << b)) != 0;
    }

    /**
     * Get a bit reference at position p for character c
     * @param c Character index [0,1,...,S[
     * @param p Position
     * @return Bit reference object
     */
    Bitref operator()(size_t c, size_t p) {
        assert(c < S);
        assert(p < N);
        size_t w = (bvSize / S) * c + (p / 64);
        size_t b = p % 64;
        return Bitref(bv[w], 1ull << b);
    }

    /**
     * Get a bit at position p for character c
     * @param c Character index [0,1,...,S[
     * @param p Position
     * @return true or false
     */
    bool at(size_t c, size_t p) const {
        return operator()(c, p);
    }

    /**
     * Create an index for the bitvector to support fast rank operations
     */
    void index() {
        // reset counts to zero
        memset((void*)counts, 0, countsSize * sizeof(size_t));

        // iterate over each character
        for (size_t c = 0; c < S; c++) {
            // the number of ones seen
            size_t numOnes = 0;
            size_t pos = (countsSize / S) * c;
            // iterate over each block  (sBlock is the start of each block)
            for (size_t sBlock = 0; sBlock < (bvSize / S);
                 sBlock += 8, pos += 2) {
                counts[pos] = numOnes;
                numOnes += __builtin_popcountll(bv[sBlock + (bvSize / S) * c]);
                for (int word = 1; word < 8; word++) {
                    counts[pos + 1] |= (numOnes - counts[pos])
                                       << 9 * (word - 1);
                    if (sBlock + word < bvSize / S) {
                        numOnes += __builtin_popcountll(
                            bv[sBlock + word + (bvSize / S) * c]);
                    }
                }
            }
            counts[(countsSize / S) * (c + 1) - 2] = numOnes;
        }
    }

    /**
     * Get the number of 1-bits within the range [0...p[ (preceding pos
     * p)
     * @param c Character index [0,1,...,S[
     * @param p Position
     */
    size_t rank(size_t c, size_t p) const {
        assert(c < S);
        assert(p < N);
        size_t w = (bvSize / S) * c + (p / 64); // word index
        const uint64_t q =
            (countsSize / S) * c + ((p / 64) / 4 & ~1); // counts index
        const int offset = (p / 64) % 8 - 1;            // offset within counts

        // add first level counts
        size_t rv = counts[q];
        // add second level counts
        rv += counts[q + 1] >>
                  (offset + (offset >> (sizeof offset * 8 - 4) & 0x8)) * 9 &
              0x1FF;
        // add popcount in final word
        return rv + __builtin_popcountll(bv[w] & ((1ULL << p % 64) - 1));
    }

    /**
     * Write the bitvector to an open filestream
     * @param ofs Open output filestream
     */
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv, bvSize * sizeof(size_t));
        ofs.write((char*)counts, countsSize * sizeof(size_t));
    }

    /**
     * Read the bitvector from an open filestream
     * @param ifs Open input filestream
     */
    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));
        allocateMem();
        ifs.read((char*)bv, bvSize * sizeof(size_t));
        ifs.read((char*)counts, countsSize * sizeof(size_t));
    }

    /**
     * Return the size of the bitvector
     * @return The size of the bitvector
     */
    size_t size() const {
        return N;
    }

    /**
     * Constructor
     * @param N Number of bits in the interleaved bitvector per
     * character
     */
    BitvecNonIntl(size_t N = 0) : N(N), bv(NULL), counts(NULL) {
        allocateMem();
    }

    /**
     * Move constructor
     * @param rhs Right hand size
     */
    BitvecNonIntl(BitvecNonIntl<S>&& rhs) : BitvecNonIntl() {
        swap(*this, rhs);
    }

    /**
     * Move assignment operator (shallow copy)
     * @param rhs Right hand size
     */
    BitvecNonIntl& operator=(BitvecNonIntl<S>&& rhs) {
        swap(*this, rhs);
        return *this;
    }

    /**
     * Deleted copy constructor and copy assignment operator
     */
    BitvecNonIntl(const BitvecNonIntl<S>&) = delete;
    BitvecNonIntl& operator=(const BitvecNonIntl<S>&) = delete;

    /**
     * Destructor
     */
    ~BitvecNonIntl() {
        free(bv);
        free(counts);
    }
};

// ============================================================================
// BIT REFERENCE CLASS
// ============================================================================

class Bitref2 {

  private:
    size_t& wordRef; // reference to a word in the bitvector
    size_t bitmask1; // bitmask of the form (01 << bitIdx)
    size_t bitmask2; // bitmask of the form (10 << bitIdx)
    size_t bitmask3; // bitmask of the form (11 << bitIdx)

  public:
    /**
     * Constructors
     * @param wordRef Reference to a word in the bitvector
     * @param bitmask Bitmask of the form (1 << bitIdx)
     */
    Bitref2(size_t& wordRef, size_t bitmask1, size_t bitmask2)
        : wordRef(wordRef), bitmask1(bitmask1), bitmask2(bitmask2) {
        bitmask3 = bitmask1 | bitmask2;
    }

    /**
     * Set bitref to particular value
     * @param val Target value
     * @return Bitref2 reference after modification
     */
    const Bitref2& operator=(int val) {
        if (val == 0) {
            wordRef &= ~bitmask3;
        } else if (val == 1) {
            wordRef &= ~bitmask3;
            wordRef |= bitmask1;
        } else if (val == 2) {
            wordRef &= ~bitmask3;
            wordRef |= bitmask2;
        } else if (val == 3) {
            wordRef |= bitmask3;
        } else {
            throw std::runtime_error(
                "New value in Bitvec2 must be between 0-3");
        }
        return *this;
    }

    /**
     * Set bitref to another bitref
     * @param br Another bitref
     * @return Bitref2 reference
     */
    const Bitref2& operator=(const Bitref2& br) {
        return this->operator=(int(br));
    }

    /**
     * Bool conversion operator
     */
    operator int() const {
        if ((wordRef & bitmask3) == 0) {
            return 0;
        } else {
            int mod = (wordRef & bitmask3) % 3;
            if (mod == 0) {
                return 3;
            } else {
                return mod;
            }
        }
    }
};

// ============================================================================
// BIT VECTOR CLASS
// ============================================================================

class Bitvec2 {

  private:
    size_t N;               // size of the bitvector
    std::vector<size_t> bv; // actual bitvector
    // std::vector<size_t> counts; // interleaved 1st and 2nd level counts

  public:
    /**
     * Get a bit at a certain position
     * @param p Position
     * @return true or false
     */
    size_t operator[](size_t p) const {
        assert(p < N);
        size_t w = 2 * p / 64;
        size_t b = 2 * p % 64;
        return (bv[w] >> b & 3ull);
    }

    /**
     * Get a bit reference at a certain position
     * @param p Position
     * @return Bit reference object
     */
    Bitref2 operator[](size_t p) {
        assert(p < N);
        size_t w = 2 * p / 64;
        size_t b = 2 * p % 64;
        return Bitref2(bv[w], 1ull << b, 2ull << b);
    }

    /**
     * Write the bitvector to an open filestream
     * @param ofs Open output filestream
     */
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv.data(), bv.size() * sizeof(size_t));
        // ofs.write((char*)counts.data(), counts.size() * sizeof(size_t));
    }

    /**
     * Read the bitvector from an open filestream
     * @param ifs Open input filestream
     */
    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));

        bv.resize(2 * (N + 63) / 64);
        ifs.read((char*)bv.data(), bv.size() * sizeof(size_t));

        // counts.resize((bv.size() + 7) / 4);
        // ifs.read((char*)counts.data(), counts.size() * sizeof(size_t));
    }

    /**
     * Return the size of the bitvector
     * @return The size of the bitvector
     */
    size_t size() const {
        return N;
    }

    /**
     * Default constructor, move constructor and move assignment operator
     */
    Bitvec2() : N(0){};
    Bitvec2(Bitvec2&& rhs) = default;
    Bitvec2& operator=(Bitvec2&& rhs) = default;

    /**
     * Deleted copy constructor and copy assignment operator
     */
    Bitvec2(const Bitvec2&) = delete;
    Bitvec2& operator=(const Bitvec2&) = delete;

    /**
     * Constructor
     * @param N Number of bits in the bitvector
     */
    Bitvec2(size_t N) : N(N), bv(2 * (N + 63) / 64, 0ull) {
    }
};

// ============================================================================
// BIT REFERENCE CLASS
// ============================================================================

class BitrefN {

  private:
    std::vector<uint8_t*> wordRef; // reference to a word in the bitvector
    size_t bitmask;                // bitmask of the form (1 << bitIdx)
    uint8_t shift;

  public:
    /**
     * Constructor
     * @param wordRef Reference to a word in the bitvector
     * @param bitmask Bitmask of the form (1 << bitIdx)
     */
    BitrefN(std::vector<uint8_t*> wordRef, size_t bitmask, uint8_t shift)
        : wordRef(wordRef), bitmask(bitmask), shift(shift) {
    }

    /**
     * Set bitref to particular value
     * @param val Target value
     * @return Bitref reference after modification
     */
    const BitrefN& operator=(size_t val) {
        *wordRef[0] &= ~(bitmask << shift);
        *wordRef[0] |= val << shift;
        for (size_t i = 1; i < wordRef.size(); i++) {
            *wordRef[i] &= ~(bitmask >> (8 - shift + (i - 1) * 8));
            *wordRef[i] |= val >> (8 - shift + (i - 1) * 8);
        }
        return *this;
    }

    /**
     * Set bitref to another bitref
     * @param br Another bitref
     * @return Bitref reference
     */
    const BitrefN& operator=(const BitrefN& br) {
        return this->operator=(int(br));
    }

    /**
     * Bool conversion operator
     */
    operator int() const {
        if (wordRef.empty()) {
            return 0;
        }
        int result = (*wordRef[0] & (bitmask << shift)) >> shift;
        for (size_t i = 1; i < wordRef.size(); i++) {
            uint8_t currentShift = (8 - shift + (i - 1) * 8);
            result += (*wordRef[i] & (bitmask >> currentShift)) << currentShift;
        }
        return result;
    }
};

class BitrefNConst {

  private:
    std::vector<const uint8_t*> wordRef; // reference to a word in the bitvector
    size_t bitmask;                      // bitmask of the form (1 << bitIdx)
    uint8_t shift;

  public:
    /**
     * Constructor
     * @param wordRef Reference to a word in the bitvector
     * @param bitmask Bitmask of the form (1 << bitIdx)
     */
    BitrefNConst(std::vector<const uint8_t*> wordRef, size_t bitmask,
                 uint8_t shift)
        : wordRef(wordRef), bitmask(bitmask), shift(shift) {
    }

    /**
     * Bool conversion operator
     */
    operator int() const {
        if (wordRef.empty()) {
            return 0;
        }
        int result = (*wordRef[0] & (bitmask << shift)) >> shift;
        for (size_t i = 1; i < wordRef.size(); i++) {
            uint8_t currentShift = (8 - shift + (i - 1) * 8);
            result += (*wordRef[i] & (bitmask >> currentShift)) << currentShift;
        }
        return result;
    }
};

// ============================================================================
// BIT VECTOR CLASS
// ============================================================================

class BitvecN {

  private:
    size_t N;    // number of elements
    uint8_t len; // length of entries
    size_t bitmask;
    std::vector<uint8_t> bv; // actual bitvector
    // std::vector<size_t> counts; // interleaved 1st and 2nd level counts

  public:
    /**
     * Get a bit reference at a certain position
     * @param p Position
     * @return Bit reference object
     */
    BitrefN operator[](size_t p) {
        assert(p < N);
        if (len == 0) {
            return BitrefN(std::vector<uint8_t*>(), bitmask, 0);
        }
        size_t w = p * len / 8;
        uint8_t b = p * len % 8;
        std::vector<uint8_t*> words = {&bv[w]};
        size_t extra = std::ceil((float)(len - (8 - b)) / (float)8);
        for (size_t i = 1; i <= extra; i++) {
            words.emplace_back(&bv[w + i]);
        }
        return BitrefN(words, bitmask, b);
    }

    /**
     * Get a bit reference at a certain position
     * @param p Position
     * @return Bit reference object
     */
    const BitrefNConst operator[](size_t p) const {
        assert(p < N);
        if (len == 0) {
            return BitrefNConst(std::vector<const uint8_t*>(), bitmask, 0);
        }
        size_t w = p * len / 8;
        uint8_t b = p * len % 8;
        std::vector<const uint8_t*> words = {&bv[w]};
        size_t extra = std::ceil((float)(len - (8 - b)) / (float)8);
        for (size_t i = 1; i <= extra; i++) {
            words.emplace_back(&bv[w + i]);
        }
        return BitrefNConst(words, bitmask, b);
    }

    // /**
    //  * Create an index for the bitvector to support fast rank operations
    //  */
    // void index() {
    //     counts = std::vector<size_t>((bv.size() + 7) / 4, 0ull);

    //     size_t countL1 = 0, countL2 = 0;
    //     for (size_t w = 0, q = 0; w < bv.size(); w++) {
    //         if (w % 8 == 0) { // store the L1 counts
    //             countL1 += countL2;
    //             counts[q] = countL1;
    //             countL2 = __builtin_popcountll(bv[w]);
    //             q += 2;
    //         } else { // store the L2 counts
    //             counts[q - 1] |= (countL2 << (((w % 8) - 1) * 9));
    //             countL2 += __builtin_popcountll(bv[w]);
    //         }
    //     }
    // }

    // /**
    //  * Get the number of 1-bits within the range [0...p[ (preceding pos p)
    //  * @param p Position
    //  */
    // size_t rank(size_t p) const {
    //     assert(p < N);
    //     size_t w = p / 64;      // word index
    //     size_t b = p % 64;      // bit offset
    //     size_t q = (w / 8) * 2; // counts index

    //     // add the first-level counts
    //     size_t rv = counts[q];

    //     // add the second-level counts
    //     int64_t t = (w % 8) - 1;
    //     rv += counts[q + 1] >> (t + (t >> 60 & 8)) * 9 & 0x1FF;

    //     // add the popcount in the final word
    //     return rv + __builtin_popcountll((bv[w] << 1) << (63 - b));
    // }

    /**
     * Write the bitvector to an open filestream
     * @param ofs Open output filestream
     */
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&len, sizeof(len));
        ofs.write((char*)bv.data(), bv.size() * sizeof(uint8_t));
        // ofs.write((char*)counts.data(), counts.size() * sizeof(size_t));
    }

    /**
     * Read the bitvector from an open filestream
     * @param ifs Open input filestream
     */
    void read(std::ifstream& ifs, size_t size) {
        N = size;

        ifs.read((char*)&len, sizeof(len));

        bv.resize((N * len + 7) / 8);
        ifs.read((char*)bv.data(), bv.size() * sizeof(uint8_t));

        bitmask = std::pow(2, len) - 1;

        // counts.resize((bv.size() + 7) / 4);
        // ifs.read((char*)counts.data(), counts.size() * sizeof(size_t));
    }

    /**
     * Return the number of elements of the bitvector
     * @return The number of elements of the bitvector
     */
    size_t nrOfElements() const {
        return N;
    }

    /**
     * Return the length of the elements of the bitvector
     * @return The length of the elements of the bitvector
     */
    uint8_t elementLength() const {
        return len;
    }

    /**
     * Return the size of the bitvector
     * @return The size of the bitvector
     */
    size_t size() const {
        return N * len;
    }

    /**
     * Default constructor, move constructor and move assignment operator
     */
    BitvecN() : N(0), len(0){};

    // // TODO do we need the code below? It gives errors
    // BitvecN(BitvecN&& rhs) = default;
    // BitvecN& operator=(BitvecN&& rhs) = default;

    // /**
    //  * Deleted copy constructor and copy assignment operator
    //  */
    // BitvecN(const BitvecN&) = delete;
    // BitvecN& operator=(const BitvecN&) = delete;

    /**
     * Constructor
     * @param N Number of bits in the bitvector
     */
    BitvecN(size_t N)
        : N(N), len(std::ceil(std::log2(N))), bitmask(std::pow(2, len) - 1),
          bv((N * len + 7) / 8, 0) {
    }
};

#endif