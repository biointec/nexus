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

#include <fstream> // I/O

// ============================================================================
// STRUCT MAPPINGPAIR
// ============================================================================

struct MappingPair {
    // identifier of the node
    int id;
    // distance to the right end of the node
    int distanceFromRightEnd;

    /**
     * @brief Construct a new Mapping Pair object
     *
     */
    MappingPair() : id(-1), distanceFromRightEnd(-1) {
    }

    /**
     * @brief Construct a new Mapping Pair object
     *
     * @param id distance to the right end of the node
     * @param distanceFromRightEnd
     */
    MappingPair(int id, int distanceFromRightEnd)
        : id(id), distanceFromRightEnd(distanceFromRightEnd) {
    }

    /**
     * @brief Writes the node to a file
     *
     * @param ofs output stream
     * @return std::ofstream& - output stream
     */
    std::ofstream& write(std::ofstream& ofs) {
        ofs.write((char*)&id, sizeof(id));
        ofs.write((char*)&distanceFromRightEnd, sizeof(distanceFromRightEnd));
        return ofs;
    }

    /**
     * @brief Reads the node from a file
     *
     * @param ifs input stream
     * @return std::ifstream& - input stream
     */
    std::ifstream& read(std::ifstream& ifs) {
        ifs.read((char*)&id, sizeof(id));
        ifs.read((char*)&distanceFromRightEnd, sizeof(distanceFromRightEnd));
        return ifs;
    }
};