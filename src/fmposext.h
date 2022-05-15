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

#include "fmocc.h"

// ============================================================================
// CLASS FMPosExt
// ============================================================================

/**
 * A single node in the bidirectional FM-index. Its depth is the depth from the
 * startmatch for a particular phase of a search
 */

template <class positionClass> class FMPosExt : public positionClass {

  private:
    // the character of this node
    char c;
    // has this particular node already reported?
    bool reported = false;

  public:
    /**
     * @brief Construct a new FMPosExt object or a node of the search tree
     *
     * @param character the character of this node
     * @param ranges the ranges over the suffix and reversed suffix array that
     * go to this node
     * @param row the row of this node in the alignment matrix = depth of this
     * node
     * @param trueDepth the true depth of the total current match
     */
    FMPosExt(char character, SARangePair ranges, length_t row,
             int trueDepth = -1)
        : positionClass(ranges, row, trueDepth), c(character), reported(false) {
    }

    /**
     * @brief Construct a new FMPosExt object or a node of the search tree
     *
     * @param character the character of this node
     * @param position the position in the FM-index or the compressed de Bruijn
     * graph on which this FMPosExt is built
     */
    FMPosExt(char character, const positionClass& position)
        : positionClass(position), c(character), reported(false) {
    }

    /**
     * @brief Construct a new FMPosExt object. Default constructor, this Node
     * will have empty ranges.
     *
     */
    FMPosExt() : positionClass(), c(char(0)) {
    }

    /**
     * @brief Set the report flag to true
     *
     */
    void report() {
        reported = true;
    }

    /**
     * @brief Report the match (with added depth) at this node
     *
     * @param occ the match will be stored here
     * @param startDepth the depth to add to the match
     * @param EDFound the found edit distance for this node
     * @param noDoubleReports false if this node is allowed to report more than
     * once, defaults to false
     * @param shift right shift of the matches, defaults to zero
     */
    void report(FMOcc<positionClass>& occ, const length_t& startDepth,
                const length_t& EDFound, const bool& noDoubleReports = false,
                length_t shift = 0) {
        if (!reported) {
            // Set the depth to the total depth of the match instead of the row
            // of this node in the alignment matrix
            this->setDepth(this->depth + startDepth);
            // Initialize the resulting FMOcc object
            occ = FMOcc<positionClass>(*this, EDFound, shift);
            // Reset the depth
            this->setDepth(this->depth - startDepth);

            // if finalPiece, report only once
            if (noDoubleReports) {
                report();
            }
        }
    }

    /**
     * @brief Get the ranges of this node
     *
     * @return const SARangePair& - the ranges of this node
     */
    const SARangePair& getRanges() const {
        return this->ranges;
    }

    /**
     * @brief Get the character of this node
     *
     * @return const char - the character of this node
     */
    const char getCharacter() const {
        return c;
    }

    /**
     * @brief Get the row of this node
     *
     * @return unsigned int - the row of this node
     */
    unsigned int getRow() const {
        return this->depth;
    }
};