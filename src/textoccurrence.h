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
#include <vector>

// ============================================================================
// CLASS TextOccurrence
// ============================================================================

class TextOccurrence {
  private:
    Range range;        // the range in the text
    int distance;       // the distance to this range (edit or hamming)
    std::string output; // the corresponding output for this occurrence (for now
                        // a custom format)

  public:
    /**
     * @brief Construct a new TextOccurrence object
     *
     * @param range, the range of this occurrence in the text
     * @param distance, the (edit or hamming) distance to the mapped read of
     * this occurrence
     */
    TextOccurrence(Range range, int distance)
        : range(range), distance(distance), output() {
    }

    /**
     * @brief Generates the output of this occurrence, for now in format:
     * startposition\twidth\tdistance, where startposition is the beginning of
     * the text occurrence, width is the length of this occurrence, distance is
     * the (edit or hamming) distance to the mapped read.
     *
     */
    void generateOutput() {
        output = std::to_string(range.getBegin()) + "\t" +
                 std::to_string(range.width()) + "\t" +
                 std::to_string(distance);
    }

    /**
     * @brief Get the range in the text
     *
     * @return const Range& - the range in the text
     */
    const Range& getRange() const {
        return range;
    }

    /**
     * @brief Get the distance to the range (edit or hamming)
     *
     * @return const int& - the distance to the range (edit or hamming)
     */
    const int& getDistance() const {
        return distance;
    }

    /**
     * @brief Get the corresponding output for this occurrence
     *
     * @return const std::string& - the corresponding output for this occurrence
     */
    const std::string& getOutput() const {
        return output;
    }

    /**
     * Operator overloading for sorting the occurrences.
     * Occurrences are first sorted on their begin position, then on their
     * distance and finally on their width
     */

    /**
     * @brief Operator overloading for sorting the occurrences. Occurrences are
     * first sorted on their begin position, then on their distance and finally
     * on their width.
     *
     * @param r Other TextOccurrence object for comparison
     * @return true if this object is smaller than the other
     * @return false otherwise
     */
    bool operator<(const TextOccurrence& r) {
        if (range.getBegin() != r.getRange().getBegin()) {
            return range.getBegin() < r.getRange().getBegin();
        } else {
            // begin is equal, better ed is smarter
            if (distance != r.getDistance()) {
                return distance < r.getDistance();
            } else {
                // shorter read is smaller...
                return range.width() < r.getRange().width();
            }
        }
    }

    /**
     * @brief Operator overloading for comparing the occurrences. Occurrences
     * are equal when their ranges and distances are equal.
     *
     * @param r Other TextOccurrence object for comparison
     * @return true if equal
     * @return false if not equal
     */
    bool operator==(const TextOccurrence& r) {
        return r.getRange() == range && r.getDistance() == distance;
    }
};

// ============================================================================
// CLASS TextOccurrenceSFI
// ============================================================================

class TextOccurrenceSFI : public TextOccurrence {

    // This is class is used to contain occurrences in the original text along
    // with their location in the graph. This class is only used during
    // strain-fixed matching.

  private:
    // The node path in the graph that corresponds to this TextOccurrence
    std::vector<uint32_t> nodepath;
    // The strain in the pan-genome to which this TextOccurrence belongs
    int strain;
    // The distance from the left end of the start node to the start of the
    // occurrence
    uint32_t distanceFromLeftEnd;

  public:
    /**
     * @brief Construct a new TextOccurrenceSFI object
     *
     * @param range the range in the text
     * @param distance the distance to this range (edit or hamming)
     * @param nodepath the node path in the graph
     * @param strain the strain in the pan-genome
     * @param distanceFromLeftEnd The distance from the left end of the start
     * node to the start of the occurrence
     */
    TextOccurrenceSFI(Range range, int distance, std::vector<uint32_t> nodepath,
                      uint32_t distanceFromLeftEnd)
        : TextOccurrence(range, distance), nodepath(nodepath), strain(-1),
          distanceFromLeftEnd(distanceFromLeftEnd) {
    }

    /**
     * @brief Get the node path in the graph
     *
     * @return const std::vector<uint32_t>& - the node path in the graph
     */
    const std::vector<uint32_t>& getNodePath() const {
        return nodepath;
    }

    /**
     * @brief Get the strain in the pan-genome
     *
     * @return const int& - the strain in the pan-genome
     */
    const int& getStrain() const {
        return strain;
    }

    /**
     * @brief Set the strain this occurrence stems from
     *
     * @param strain
     */
    void setStrain(const int& strain) {
        this->strain = strain;
    }

    /**
     * @brief Get the distance from the left end of the first node in the
     * pan-genome
     *
     * @return const uint32_t& - the distance from the left end of the first
     * node in the pan-genome
     */
    const uint32_t& getDistanceFromLeftEnd() const {
        return distanceFromLeftEnd;
    }
};