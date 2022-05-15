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

#pragma once

#include "range.h"

// ============================================================================
// CLASS SARANGEPAIR
// ============================================================================

/**
 * A pair of ranges. The first range is range over the suffix array of the text.
 * The second range is the corresponding range over the suffix array of the
 * reversed text
 */

class SARangePair {
  private:
    Range rangeSA;    // the range over the suffix array
    Range rangeSARev; // the range over the suffix array of the reversed text

  public:
    /**
     * @brief Construct a new SARangePair object. Default constructor, creates
     * two empty ranges.
     *
     */
    SARangePair() : rangeSA(Range()), rangeSARev(Range()) {
    }

    /**
     * @brief Construct a new SARangePair object.
     *
     * @param rangeSA the range over the suffix array
     * @param rangeSARev the range over the suffix array of the reversed text
     */
    SARangePair(Range rangeSA, Range rangeSARev)
        : rangeSA(rangeSA), rangeSARev(rangeSARev) {
    }

    /**
     * @brief Get the range over the suffix array
     *
     * @return const Range& - the range over the suffix array
     */
    const Range& getRangeSA() const {
        return rangeSA;
    }

    /**
     * @brief Get the range over the suffix array of the reversed text
     *
     * @return const Range& - the range over the suffix array of the reversed
     * text
     */
    const Range& getRangeSARev() const {
        return rangeSARev;
    }

    /**
     * @brief Check if the ranges are empty.
     *
     * @return true if the ranges are empty
     * @return false otherwise
     */
    bool empty() const {
        return rangeSA.empty();
    }

    /**
     * @brief Get the width of the ranges.
     *
     * @return length_t - the width of the ranges
     */
    length_t width() const {
        return rangeSA.width();
    }

    /**
     * @brief Operator overloading to check if two SARangePair objects are
     * equal.
     *
     * @param o Other SARangePair object for comparison
     * @return true if equal
     * @return false if not equal
     */
    bool operator==(const SARangePair& o) const {
        // only the first range matters as the ranges imply each other
        return o.getRangeSA() == rangeSA;
    }
};