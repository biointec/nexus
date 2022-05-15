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

#include "bitvec.h"

#include <stdint.h>

// ============================================================================
// (TYPE) DEFINITIONS AND PROTOTYPES
// ============================================================================

typedef uint64_t length_t;

// ============================================================================
// STRUCT NODE
// ============================================================================

struct Node {
    // length of the substring corresponding to the node
    uint32_t len;
    // multiplicity of the node
    uint32_t multiplicity;
    // left boundary of the suffix interval of the substring corresponding to
    // the node
    length_t left_kmer_forward;
    // left boundary of the suffix interval of the k-length suffix of the
    // substring corresponding to the node
    length_t right_kmer_forward;
    // left boundary of the suffix interval in the reverse SA of the substring
    // corresponding to the node
    length_t right_kmer_reverse;
    // left boundary of the suffix interval in the reverse SA of the k-length
    // suffix of the substring corresponding to the node
    length_t left_kmer_reverse;
    // Mapping of regular and reverse ranks of edges passing through the node.
    // Specifically, regular ranks are mapped to revers ranks.
    BitvecN edgeMapping;
    // for visualization
    bool visited = false;

    /**
     * @brief Construct a new Node object
     *
     */
    Node()
        : len(0), multiplicity(0), left_kmer_forward(0), right_kmer_forward(0),
          right_kmer_reverse(0), left_kmer_reverse(0), edgeMapping(0) {
    }

    /**
     * @brief Construct a new Node object
     *
     * @param len, length of the substring corresponding to the node
     * @param multiplicity, multiplicity of the node
     * @param left_kmer_forward, left boundary of the suffix interval of the
     * substring corresponding to the node
     * @param right_kmer_forward, left boundary of the suffix interval of the
     * k-length suffix of the substring corresponding to the node
     * @param right_kmer_reverse, left boundary of the suffix interval in the
     * reverse SA of the substring corresponding to the node
     * @param left_kmer_reverse, left boundary of the suffix interval in the
     * reverse SA of the k-length suffix of the substring corresponding to the
     * node
     */
    Node(int len, int multiplicity, length_t left_kmer_forward,
         length_t right_kmer_forward, length_t right_kmer_reverse,
         length_t left_kmer_reverse)
        : len(len), multiplicity(multiplicity),
          left_kmer_forward(left_kmer_forward),
          right_kmer_forward(right_kmer_forward),
          right_kmer_reverse(right_kmer_reverse),
          left_kmer_reverse(left_kmer_reverse), edgeMapping(multiplicity) {
    }

    /**
     * @brief Writes the node to a file
     *
     * @param ofs output strean
     * @return std::ofstream& - output stream
     */
    std::ofstream& write(std::ofstream& ofs) {
        ofs.write((char*)&len, sizeof(len));
        ofs.write((char*)&multiplicity, sizeof(multiplicity));
        ofs.write((char*)&left_kmer_forward, sizeof(left_kmer_forward));
        ofs.write((char*)&right_kmer_forward, sizeof(right_kmer_forward));
        ofs.write((char*)&right_kmer_reverse, sizeof(right_kmer_reverse));
        ofs.write((char*)&left_kmer_reverse, sizeof(left_kmer_reverse));
        edgeMapping.write(ofs);
        return ofs;
    }

    /**
     * @brief Reads the node from a file
     *
     * @param ifs input stream
     * @return std::ifstream& - input stream
     */
    std::ifstream& read(std::ifstream& ifs) {
        ifs.read((char*)&len, sizeof(len));
        ifs.read((char*)&multiplicity, sizeof(multiplicity));
        ifs.read((char*)&left_kmer_forward, sizeof(left_kmer_forward));
        ifs.read((char*)&right_kmer_forward, sizeof(right_kmer_forward));
        ifs.read((char*)&right_kmer_reverse, sizeof(right_kmer_reverse));
        ifs.read((char*)&left_kmer_reverse, sizeof(left_kmer_reverse));
        edgeMapping.read(ifs, multiplicity);
        visited = false;
        return ifs;
    }
};

struct VisualizationNode : Node {
    // Stores if this node is part of the original node path
    bool part_of_path;
    // The string corresponding to the node
    std::string omega;
    // The short version of the string corresponding to the node
    std::string omega_short;

    /**
     * @brief Construct a new Visualization Node from an existing node and its
     * additional attributes
     *
     * @param node The original node in the graph
     * @param part_of_path Stores if this node is part of the original node path
     * @param omega The string corresponding to the node
     * @param omega_short The short version of the string corresponding to the
     * node
     */
    VisualizationNode(Node& node, bool& part_of_path, std::string& omega,
                      std::string& omega_short)
        : Node(node.len, node.multiplicity, node.left_kmer_forward,
               node.right_kmer_forward, node.right_kmer_reverse,
               node.left_kmer_reverse),
          part_of_path(part_of_path), omega(omega), omega_short(omega_short) {
    }
};