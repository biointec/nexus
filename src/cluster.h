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

#include "fmposext.h"

// ============================================================================
// CLASS CLUSTER
// ============================================================================

template <class positionClass> class Cluster {
  private:
    // the edit distances of this cluster
    std::vector<uint> eds;
    // the nodes of this cluster
    std::vector<FMPosExt<positionClass>> nodes;

    // the lastCell of the cluster that was filled in
    uint lastCell;
    // the maxEd for this cluster
    uint maxED;
    // the startdepth for this cluster (= depth of match before matrix of this
    // cluster)
    uint startDepth;

    // the right shift of the occurrences in the text
    uint shift;

  public:
    /**
     * @brief Construct a new Cluster object
     *
     * @param size the size of the cluster
     * @param maxED the maximal allowed edit distance
     * @param startDepth the depth before this cluster
     * @param shift the right shift of the occurrences in the text
     */
    Cluster(int size, int maxED, int startDepth, int shift)
        : eds(size, maxED + 1), nodes(size), lastCell(-1), maxED(maxED),
          startDepth(startDepth), shift(shift) {
    }

    /**
     * @brief Set the ed and node at index idx to ed and node. Also updates
     * lastCell to be idx
     *
     * @param idx the idx to change
     * @param node the node to set at index idx
     * @param ed the ed to set at index idx
     */
    void setValue(int idx, const FMPosExt<positionClass>& node, const int& ed) {
        eds[idx] = ed;
        nodes[idx] = node;
        lastCell = idx;
    }

    /**
     * @brief Return the size of this cluster
     *
     * @return const unsigned int - the size of this cluster
     */
    const unsigned int size() const {
        return eds.size();
    }

    /**
     * @brief Fills in a vector with all nodes in the cluster that are a centre
     * and under the maximal allowed distance, be aware that if there are
     * multiple centers in the cluster it is very likely that only one of them
     * will be redundant, but the others might eliminate another occurrence.
     * This function is only used for strain-fixed mapping.
     *
     * @param centers - the vector of center nodes
     */
    void reportCentersAtEnd(std::vector<FMOcc<FMPos>>& centers,
                            std::vector<uint32_t>& leftNodes,
                            std::vector<uint32_t>& rightNodes) {

        for (uint i = 0; i <= lastCell; i++) {
            if (eds[i] <= maxED && (i == 0 || eds[i] <= eds[i - 1]) &&
                (i == lastCell || eds[i] <= eds[i + 1])) {
                FMOcc<FMPos> m;
                nodes[i].report(m, startDepth, eds[i], true, shift);
                centers.emplace_back(m);
            }
        }
    }

    /**
     * @brief Fills in a vector with all nodes in the cluster that are a centre
     * and under the maximal allowed distance. Centers are searched within all
     * nodes that have the same node path. Be aware that if there are multiple
     * centers in the cluster it is very likely that only one of them will be
     * redundant, but the others might eliminate another occurrence. This
     * function is only used for strain-free mapping.
     *
     * @param centers - the vector of center nodes
     */
    void reportCentersAtEnd(std::vector<FMOcc<FMPosSFR>>& centers,
                            std::vector<uint32_t>& leftNodes,
                            std::vector<uint32_t>& rightNodes) {
        // Find all start positions of sets of nodes that have the same node
        // path.
        std::vector<uint> startpositions = {0};

        auto& hf_nodes = static_cast<std::vector<FMPosExt<FMPosSFR>>&>(nodes);

        const auto& firstPos = hf_nodes[0];

        uint32_t prevSize = firstPos.nodePathSize();

        for (uint i = 1; i <= lastCell; i++) {
            const auto& pos = hf_nodes[i];

            uint32_t newSize = pos.nodePathSize();

            if (newSize > prevSize) {
                startpositions.emplace_back(i);
                prevSize = newSize;
            }
        }
        // Also add the end of the last set of nodes
        startpositions.emplace_back(lastCell + 1);

        // Find the cluster centres for every set of nodes with the same node
        // path separately.
        for (uint j = 0; j < startpositions.size() - 1; j++) {
            for (uint i = startpositions[j]; i < startpositions[j + 1]; i++) {
                if (eds[i] <= maxED &&
                    (i == startpositions[j] || eds[i] <= eds[i - 1]) &&
                    (i == startpositions[j + 1] - 1 || eds[i] <= eds[i + 1])) {
                    FMOcc<FMPosSFR> m;

                    auto& pos = hf_nodes[i];

                    // find nodePath at this particular node
                    // start by making a copy of left
                    auto path = leftNodes;

                    // constrict left to  appropriate number of nodes
                    path.resize(pos.getNumberOfNodesLeft());

                    // reverse left
                    std::reverse(path.begin(), path.end());
                    // add right to left
                    path.insert(path.end(), rightNodes.begin(),
                                rightNodes.begin() +
                                    pos.getNumberOfNodesRight());

                    pos.setNodePath(path);
                    pos.report(m, startDepth, eds[i], true, shift);
                    centers.emplace_back(m);
                }
            }
        }
    }

    /**
     * @brief Returns the approximate match that corresponds to the ranges of
     * the deepest global minimum of this cluster, but with the depth of the
     * highest global minimum of this cluster. If the direction is backward, a
     * shift will be set such that the occurrence in the text will be as short
     * as possible.
     *
     * @param dir The current direction
     * @return FMOcc<positionClass> - the approximate match that corresponds to
     * the ranges of the deepest global minimum of this cluster, but with the
     * depth of the highest global minimum of this cluster
     */
    FMOcc<positionClass> reportDeepestMinimum(Direction dir) {
        uint minED = maxED + 1;
        int highestBestIdx = -1;
        int deepestBestIdx = -1;

        for (uint i = 0; i <= lastCell; i++) {
            if (eds[i] < minED) {
                minED = eds[i];
                highestBestIdx = i;
                deepestBestIdx = i;
            }
            if (eds[i] == minED) {
                deepestBestIdx = i;
            }
        }
        FMOcc<positionClass> m;
        if (minED <= maxED) {
            nodes[deepestBestIdx].report(
                m, startDepth - (deepestBestIdx - highestBestIdx), minED, true,
                ((dir == BACKWARD) ? (deepestBestIdx - highestBestIdx) : 0) +
                    shift);
        }
        return m;
    }

    /**
     * @brief Report the matches for which the search procedure must be
     * continued when the end of the pattern is reached (left or right), but the
     * search is yet to be continued on the other side. This function is only
     * used during strain-fixed matching.
     *
     * @param dir The current direction
     * @param matches A vector of approximate matches for which the search needs
     * to be continued in the other direction. This vector will be filled in. It
     * will contain only one match.
     */
    void reportNonFinalEdgePiece(Direction dir,
                                 std::vector<FMOcc<FMPos>>& matches) {
        matches.emplace_back(reportDeepestMinimum(dir));
    }

    /**
     * @brief Report the matches for which the search procedure must be
     * continued when the end of the pattern is reached (left or right), but the
     * search is yet to be continued on the other side. This function is only
     * used during strain-free matching.
     *
     * @param dir The current direction
     * @param matches A vector of approximate matches for which the search needs
     * to be continued in the other direction. This vector will be filled in.
     */

    void reportNonFinalEdgePiece(Direction dir,
                                 std::vector<FMOcc<FMPosSFR>>& matches) {
        uint i = 0;

        auto& hf_nodes = dynamic_cast<std::vector<FMPosExt<FMPosSFR>>&>(nodes);

        while (i <= lastCell) {

            // Nodepaths of children positions either have one node more than or
            // are equal to their parent's nodepath Knowing if the nodepath is
            // different boils down to checking if the size of the nodepath is
            // different
            uint32_t prevSize = hf_nodes[0].nodePathSize(), newSize = prevSize;

            uint minED = maxED + 1;
            int highestBestIdx = -1;
            // If the node representation has not been changed yet, report all
            // (we cannot know yet which occurrences have different node paths)
            while (newSize == 0) {
                if (eds[i] <= maxED) {
                    FMOcc<positionClass> m;
                    nodes[i].report(m, startDepth, eds[i], true, shift);
                    matches.push_back(m);
                }
                i++;
                if (i <= lastCell) {
                    newSize = hf_nodes[i].nodePathSize();
                } else {
                    return;
                }
            }
            // Get one match for every set of nodes with the same node path
            while (prevSize == newSize && i <= lastCell) {
                if (eds[i] < minED) {
                    minED = eds[i];
                    highestBestIdx = i;
                }
                i++;
                if (i <= lastCell) {
                    newSize = hf_nodes[i].nodePathSize();
                }
            }
            if (minED <= maxED) {
                FMOcc<positionClass> m;
                nodes[highestBestIdx].report(m, startDepth, eds[highestBestIdx],
                                             true, shift);
                matches.emplace_back(m);
            }
        }
    }

    /**
     * @brief This method returns a match that corresponds to the highest
     * cluster centre. Its descendants and the corresponding initialization eds
     * are updated. Eds of descendants that are part of a cluster centre which
     * is lower than the lowerbound will be updated in the initEds vector
     *
     * @param lowerBound the lowerbound for this iteration
     * @param desc the descendants of the highest cluster centre, these will be
     * inserted during the method
     * @param initEds the initialization eds for the next iteration, these
     * correspond to the eds of the highest centre and its descendants, where
     * eds part of a cluster of which the centre is below the lowerbound are
     * updated. These values will be inserted during the method
     * @return FMOcc<positionClass> - The occurrence corresponding to the upper
     * cluster centre which has a valid distance
     */
    FMOcc<positionClass>
    getClusterCentra(uint lowerBound,
                     std::vector<FMPosExt<positionClass>>& desc,
                     std::vector<uint>& initEds) {
        desc.reserve(eds.size());
        initEds.reserve(eds.size());
        FMOcc<positionClass> m;
        for (uint i = 0; i <= lastCell; i++) {
            if (eds[i] > maxED || eds[i] < lowerBound) {
                continue;
            }
            bool betterThanParent = (i == 0) || eds[i] <= eds[i - 1];
            bool betterThanChild = (i == lastCell) || eds[i] <= eds[i + 1];

            if (betterThanParent && betterThanChild) {
                // this is a valid centre
                nodes[i].report(m, startDepth, eds[i], false, shift);

                // get all the descendants
                initEds.emplace_back(eds[i]);
                for (uint j = i + 1; j <= lastCell; j++) {
                    desc.emplace_back(nodes[j]);
                    initEds.emplace_back(eds[j]);
                }

                // replace the clusters under the lowerbound
                for (unsigned int k = 1; k < initEds.size(); k++) {
                    if (initEds[k] < lowerBound &&
                        initEds[k] <= initEds[k - 1] &&
                        (k == initEds.size() - 1 ||
                         initEds[k] <= initEds[k + 1])) {
                        // k is a centre under the lowerbound

                        unsigned int highestPoint = 0;
                        unsigned int lowestPoint = initEds.size() - 1;
                        // find highest point of this cluster
                        for (unsigned int l = k; l-- > 0;) {
                            if (initEds[l] != initEds[l + 1] + 1) {
                                highestPoint = l + 1;
                                break;
                            }
                        }
                        // find lowest point of this cluster
                        for (unsigned int l = k + 1; l < initEds.size(); l++) {
                            if (initEds[l] != initEds[l - 1] + 1) {
                                lowestPoint = l - 1;
                                break;
                            }
                        }

                        // highest and lowest cannot span entire initEds.size(),
                        // otherwise there would not be a valid cluster centre
                        // above the lowerbound
                        if (highestPoint != 0 &&
                            lowestPoint != initEds.size() - 1) {
                            // Make /\ with ed values of this cluster
                            // do iE[hp] = ie[hp - 1] + 1 and iE[lp] = iE[lp +
                            // 1] +1 until entire cluster has been replaced
                            length_t lC = lowestPoint;
                            length_t hC = highestPoint;
                            bool highest = true;
                            // do not go over maxED + 1, to ensure continuity at
                            // the other end
                            while (lC > hC) {
                                if (highest) {
                                    initEds[hC] = std::min(maxED + 1,
                                                           initEds[hC - 1] + 1);
                                    hC++;
                                } else {
                                    initEds[lC] = std::min(maxED + 1,
                                                           initEds[lC + 1] + 1);
                                    lC--;
                                }
                                highest = !highest;
                            }
                            if (lC == hC) {
                                // change middle element of cluster
                                initEds[lC] = std::min(initEds[lC + 1] + 1,
                                                       initEds[lC - 1] + 1);
                            }

                        } else if (highestPoint == 0 &&
                                   lowestPoint != initEds.size() - 1) {
                            // monotonous rise from lowestPoint to highestPoint
                            for (unsigned int l = lowestPoint; l-- > 0;) {
                                initEds[l] = initEds[l + 1] + 1;
                            }
                        } else if (highestPoint != 0 &&
                                   lowestPoint == initEds.size() - 1) {
                            // monotonous rise from highestPoint to lowestPoint
                            for (unsigned int l = highestPoint;
                                 l < initEds.size(); l++) {
                                initEds[l] = initEds[l - 1] + 1;
                            }
                        }
                    }
                }
                // stop searching
                break;
            }
        }

        return m;
    }
};