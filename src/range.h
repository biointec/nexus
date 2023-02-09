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

#include <ostream>

// ============================================================================
// (TYPE) DEFINITIONS AND PROTOTYPES
// ============================================================================

typedef uint64_t length_t;

// ============================================================================
// CLASS RANGE
// ============================================================================

class Range {
  private:
    length_t begin; // beginning of the range
    length_t end;   // end of the range (non-inclusive)

  public:
    /**
     * @brief Construct a new Range object
     *
     * @param b, the beginning of the range
     * @param e, the end of the range (non-inclusive)
     */
    Range(length_t b, length_t e) : begin(b), end(e) {
    }

    /**
     * @brief Construct a new Range object. Default constructor, initializes an
     * empty range.
     *
     */
    Range() : begin(0), end(0) {
    }

    /**
     * @brief Get the beginning of the range
     *
     * @return length_t - the beginning of the range
     */
    length_t getBegin() const {
        return begin;
    }

    /**
     * @brief Get the end of the range (non-inclusive)
     *
     * @return length_t - the end of the range (non-inclusive)
     */
    length_t getEnd() const {
        return end;
    }

    /**
     * @brief Check if this range is empty
     *
     * @return true if the range is empty
     * @return false otherwise
     */
    bool empty() const {
        return end <= begin;
    }

    /**
     * @brief Get the width of the range (end - begin)
     *
     * @return length_t - the width of this range
     */
    length_t width() const {
        return (empty()) ? 0 : end - begin;
    }

    /**
     * @brief Operator overloading, two ranges are equal if their begin and end
     * field are equal
     *
     * @param o Other range object for comparison
     * @return true if equal
     * @return false if not equal
     */
    bool operator==(const Range& o) const {
        return o.getBegin() == begin && o.getEnd() == end;
    }

    /**
     * @brief Check if the range contains a certain index
     *
     * @param position index to be checked
     * @return true if the index lies within the range
     * @return false otherwise
     */
    bool contains(const length_t& position) const {
        return begin <= position && position < end;
    }

    /**
     * @brief Operator overloading. Outputs the range as [begin, end) to the
     * outputstream
     *
     * @param os, the output stream
     * @param r, the range to print
     * @return std::ostream& - the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Range& r);
};

/**
 * @brief Operator overloading. Outputs the range as [begin, end) ot the
 * outputstream
 *
 * @param output, the output stream
 * @param r, the range to print
 * @return std::ostream& - the output stream
 */
std::ostream& operator<<(std::ostream& output, const Range& r);