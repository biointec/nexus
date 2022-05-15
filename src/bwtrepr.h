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

#ifndef BWTREPR_H
#define BWTREPR_H

#include <vector>

#include "alphabet.h"
#include "bitvec.h"
#include "selectinterface.h"

// ============================================================================
// CLASS BWT REPRESENTATION (supports occ(c,k) and cumOcc(c,k) in O(1) time)
// ============================================================================

template <size_t S> // S is the size of the alphabet (including '$')
class BWTRepr {     // e.g. S = 5 for DNA (A,C,G,T + $)

  private:
    // The '$' character (cIdx == 0) is not encoded in the bitvector.
    // Hence, we use only S-1 bitvectors.

    BitvecIntl<S - 1> bv; // bitvector representation of the BWT
    size_t dollarPos;     // position of the dollar sign

  public:
    /**
     * Default constructor
     */
    BWTRepr() {
    }

    /**
     * Constructor
     * @param sigma Alphabet
     * @param BWT Burrows-Wheeler transformation
     */
    BWTRepr(const Alphabet<S>& sigma, const std::string& BWT)
        : bv(BWT.size() + 1), dollarPos(BWT.size()) {
        // The $-character (cIdx == 0) is not encoded in the bitvector.
        // Hence, use index cIdx-1 in the bitvector.

        for (size_t i = 0; i < BWT.size(); i++) {
            if (BWT[i] == '$') {
                dollarPos = i;
                continue;
            }

            for (size_t cIdx = sigma.c2i(BWT[i]); cIdx < S; cIdx++)
                bv(cIdx - 1, i) = true;
        }

        bv.index();
    }

    /**
     * Get occurrence count of character c in the range BWT[0...k[
     * @param cIdx Character index
     * @param k index
     * @return occ(c, k)
     */
    size_t occ(int cIdx, size_t k) const {
        // The $-character (cIdx == 0) is not encoded in the bitvector.
        // Hence, use index cIdx-1 in the bitvector.

        if (cIdx == 0) // special case for $-character
            return (k <= dollarPos) ? 0 : 1;

        return (cIdx == 1) ? bv.rank(cIdx - 1, k)
                           : bv.rank(cIdx - 1, k) - bv.rank(cIdx - 2, k);
    }

    /**
     * Get cumulative occurrence count of characters SMALLER than c
     * in the range BWT[0...k[
     * @param cIdx Character index
     * @param k index
     * @return cumOcc(c, k)
     */
    size_t cumOcc(int cIdx, size_t k) const {
        // The $-character (cIdx == 0) is not encoded in the bitvector.
        // Hence, use index cIdx-1 in the bitvector.

        if (cIdx == 0) // special case for $-character
            return 0;

        return (cIdx == 1) ? ((k <= dollarPos) ? 0 : 1)
                           : bv.rank(cIdx - 2, k) + ((k <= dollarPos) ? 0 : 1);
    }

    /**
     * Write table to disk
     * @param filename File name
     */
    void write(const std::string& filename) {
        std::ofstream ofs(filename);
        if (!ofs)
            throw std::runtime_error("Cannot open file: " + filename);

        ofs.write((char*)&dollarPos, sizeof(dollarPos));
        bv.write(ofs);
    }

    /**
     * Load table from disk
     * @param filename File name
     */
    bool read(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs)
            return false;

        ifs.read((char*)&dollarPos, sizeof(dollarPos));
        bv.read(ifs);

        return true;
    }
};

template <size_t S>            // S is the size of the alphabet (including '$')
class BWTReprSelectSupported { // e.g. S = 5 for DNA (A,C,G,T + $)

  private:
    std::array<SelectInterface, S - 1> interfaces;
    size_t dollarPos;
    SelectOption option;

  public:
    /**
     * Constructor
     * @param sigma Alphabet
     * @param BWT Burrows-Wheeler transformation
     */
    BWTReprSelectSupported(const Alphabet<S>& sigma, const std::string& BWT,
                           const SelectOption& option)
        : dollarPos(BWT.size()), option(option) {
        // The $-character (cIdx == 0) is not encoded in the bitvector.
        // Hence, use index cIdx-1 in the bitvector.

        for (auto& i : interfaces) {
            i.setN(BWT.size());
            i.setOption(option);
        }

        for (size_t i = 0; i < BWT.size(); i++) {
            if (BWT[i] == '$') {
                dollarPos = i;
                continue;
            }

            // no prefix occurrence
            interfaces[sigma.c2i(BWT[i]) - 1][i] = true;
        }

        // index the interfaces
        for (auto& i : interfaces) {
            i.indexInterface();
        }
    }

    /**
     * Default constructor
     */
    BWTReprSelectSupported() {
    }

    /**
     * Deleted copy constructor and copy assignment operator
     */
    BWTReprSelectSupported(const BWTReprSelectSupported<S>&) = delete;
    BWTReprSelectSupported&
    operator=(const BWTReprSelectSupported<S>&) = delete;

    /**
     * Compute the index of the k-th occurrence of the character with index idx
     * in the BWT.
     * @param idx idx is the index in the count array of the character we are
     * interested in
     * @param k the number of occurrences
     * @return the index in the BWT
     */
    size_t select(size_t idx, size_t k) {
        if (idx == 0) {
            return dollarPos;
        } else {
            return interfaces[idx - 1].select(k - 1);
        }
    }

    /**
     * Write table to disk
     * @param filename File name
     */
    void write(const std::string& filename) {
        std::ofstream ofs(filename);
        if (!ofs)
            throw std::runtime_error("Cannot open file: " + filename);

        ofs.write((char*)&option, sizeof(option));

        ofs.write((char*)&dollarPos, sizeof(dollarPos));
        for (const auto& i : interfaces) {
            i.write(ofs);
        }
    }

    /**
     * Load table from disk
     * @param filename File name
     */
    bool read(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs)
            return false;

        ifs.read((char*)&option, sizeof(option));
        for (auto& i : interfaces) {
            i.setOption(option);
        }
        ifs.read((char*)&dollarPos, sizeof(dollarPos));
        for (auto& i : interfaces) {
            i.read(ifs);
        }

        return true;
    }
};

#endif