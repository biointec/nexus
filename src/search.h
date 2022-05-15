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

#include "substring.h"

#include <algorithm> // min_element
#include <cassert>
#include <vector>

// ============================================================================
// CLASS SEARCH
// ============================================================================
class Search {
  private:
    // the vector with lowerbounds
    std::vector<int> lowerBounds;
    // the vector with upperbounds
    std::vector<int> upperBounds;
    // the vector with the order of the parts
    std::vector<int> order;
    // the directions of each phase
    std::vector<Direction> directions;
    // has the direction switched for each phase
    std::vector<bool> directionSwitch;

    /**
     * @brief Construct a new Search object
     *
     * @param order the vector with the order of the parts
     * @param lowerBounds the vector with lowerbounds
     * @param upperBounds the vector with upperbounds
     * @param directions the directions of each phase
     * @param dSwitch has the direction switched for each phase
     */
    Search(std::vector<int>& order, std::vector<int>& lowerBounds,
           std::vector<int>& upperBounds, std::vector<Direction>& directions,
           std::vector<bool>& dSwitch)
        : lowerBounds(lowerBounds), upperBounds(upperBounds), order(order),
          directions(directions), directionSwitch(dSwitch) {
    }

  public:
    /**
     * @brief Static function to construct a search. The directions and switches
     * of the search are calculated
     *
     * @param order the order of the search
     * @param lowerBounds the lowerbounds of the search
     * @param upperBounds the upperbounds of the search
     * @return Search - the resulting search object
     */
    static Search makeSearch(std::vector<int> order,
                             std::vector<int> lowerBounds,
                             std::vector<int> upperBounds) {
        // check correctness of sizes
        if (order.size() != lowerBounds.size() ||
            order.size() != upperBounds.size()) {
            throw std::runtime_error("Could not create search, the sizes of "
                                     "all vectors are not equal");
        }

        // compute the directions
        std::vector<Direction> directions;
        directions.reserve(order.size());
        directions.push_back((order[1] > order[0]) ? FORWARD : BACKWARD);

        for (length_t i = 1; i < order.size(); i++) {
            Direction d = (order[i] > order[i - 1]) ? FORWARD : BACKWARD;
            directions.push_back(d);
        }

        // compute the directionswitches
        std::vector<bool> directionSwitch;
        directionSwitch.reserve(order.size());
        // first partition is not a switch
        directionSwitch.push_back(false);

        // second partition is never a switch
        // TODO check in case of Kianfar
        directionSwitch.push_back(false);

        // add the other partitions
        for (length_t i = 2; i < directions.size(); i++) {
            directionSwitch.push_back(directions[i] != directions[i - 1]);
        }

        return Search(order, lowerBounds, upperBounds, directions,
                      directionSwitch);
    }

    /**
     * @brief Set the directions of the parts to the directions of the search
     *
     * @param parts the parts to set the direction of
     */
    void setDirectionsInParts(std::vector<Substring>& parts) const {
        // set the directions for the parts
        for (length_t i = 0; i < order.size(); i++) {
            parts[order[i]].setDirection(directions[i]);
        }
    }

    /**
     * @brief Get the lowerbound for the idx'th part
     *
     * @param idx index for the part vector
     * @return int - the lowerbound for the idx'th part
     */
    int getLowerBound(int idx) const {
        assert(idx < (int)lowerBounds.size());
        return lowerBounds[idx];
    }

    /**
     * @brief Get the upperbound for the idx'th part
     *
     * @param idx index for the part vector
     * @return int - the lowerbound for the idx'th part
     */
    int getUpperBound(int idx) const {
        assert(idx < (int)upperBounds.size());
        return upperBounds[idx];
    }

    /**
     * @brief Get the idx'th part
     *
     * @param idx index for the part vector
     * @return int - the idx'th part
     */
    int getPart(int idx) const {
        assert(idx < (int)order.size());
        return order[idx];
    }

    /**
     * @brief Get the direction for the idx'th part
     *
     * @param idx index for the part vector
     * @return Direction - the direction for the idx'th part
     */
    Direction getDirection(int idx) const {
        assert(idx < (int)directions.size());
        return directions[idx];
    }

    /**
     * @brief Check if the direction switches at the idx'th part
     *
     * @param idx index for the part vector
     * @return true if the direction switches at the idx'th part
     * @return false otherwise
     */
    bool getDirectionSwitch(int idx) const {
        assert(idx < (int)directionSwitch.size());
        return directionSwitch[idx];
    }

    /**
     * @brief Get the number of parts in this search
     *
     * @return int - the number of parts in this search
     */
    int getNumParts() const {
        return order.size();
    }

    /**
     * @brief Check if the idx'th part is the first or last part of the pattern
     *
     * @param idx index for the part vector
     * @return true if the idx'th part is the first or last part of the pattern
     * @return false otherwise
     */
    bool isEdge(int idx) const {
        return order[idx] == 0 || order[idx] == getNumParts() - 1;
    }

    /**
     * @brief Check if the idx'th part is the final part of the search
     *
     * @param idx index for the part vector
     * @return true if the idx'th part is the final part of the search
     * @return false otherwise
     */
    bool isEnd(int idx) const {
        return idx == (int)order.size() - 1;
    }

    /**
     * @brief Check if the connectivity property is satisfied
     *
     * @return true if the connectivity property is satisfied
     * @return false otherwise
     */
    bool connectivitySatisfied() const {
        int highestSeen = order[0];
        int lowestSeen = order[0];
        for (length_t i = 1; i < order.size(); i++) {
            if (order[i] == highestSeen + 1) {
                highestSeen++;
            } else if (order[i] == lowestSeen - 1) {
                lowestSeen--;
            } else {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Check that the upper and lower bounds do not decrease
     *
     * @return true if the upper and lower bounds do not decrease
     * @return false otherwise
     */
    bool noDecreasingInBounds() const {
        for (length_t i = 1; i < order.size(); i++) {
            if (lowerBounds[i] < lowerBounds[i - 1]) {
                return false;
            }
            if (upperBounds[i] < upperBounds[i - 1]) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Check that the order vector starts with 0
     *
     * @return true if the order vector starts with 0
     * @return false otherwise
     */
    bool zeroBased() const {
        return *std::min_element(order.begin(), order.end()) == 0;
    }
};

/**
 * @brief Operator overloading. Outputs the search to the outputstream.
 *
 * @param os the output stream
 * @param obj the search to print
 * @return std::ostream& - the output stream
 */
std::ostream& operator<<(std::ostream& os, const Search& obj);