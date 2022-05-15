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

#include "searchstrategy.h"

class StrainFreeMapper;

// Pointer to a partition function
using PartitionPtrSFR = void (StrainFreeMapper::*)(
    const std::string&, std::vector<Substring>&, const int& numParts,
    const int& maxScore, std::vector<FMPosSFR>& exactMatchPositions,
    std::vector<std::vector<uint32_t>>& nodePaths) const;

class StrainFreeMapper {
  private:
    SearchStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>* strategy;

    PartitionPtrSFR partitionPtrSFR; // pointer to the partition method

  public:
    /**
     * @brief Construct a new strain-free mapper object
     *
     * @param strategy the search strategy that will be used
     */
    StrainFreeMapper(
        SearchStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>* strategy)
        : strategy(strategy) {

        PartitionStrategy p = strategy->getPartitioningStrategyInt();
        // set the partition strategy
        switch (p) {
        case UNIFORM:
            partitionPtrSFR = &StrainFreeMapper::partitionUniformSFR;
            break;
        case DYNAMIC:
            partitionPtrSFR = &StrainFreeMapper::partitionDynamicSFR;
            break;
        case STATIC:
            partitionPtrSFR = &StrainFreeMapper::partitionOptimalStaticSFR;
            break;
        default:
            break;
        }
    }

    /**
     * @brief Approximately match the pattern to the compressed de Bruijn graph
     * using a maximum allowed edit distance of maxED. This is done in a
     * strain-free way.
     *
     * @param pattern The pattern to be matched
     * @param maxED The maximum allowed edit distance
     * @return std::vector<FMOccSFR> - a vector of resulting occurrences
     */
    std::vector<FMOccSFR> matchApproxSFR(const std::string& pattern,
                                         length_t maxED) const;

    /**
     * @brief Executes the search recursively. If U[0] != 1, then the search
     * will start at pi[0], else the search will start with idx i and U[i]!=0
     * and U[j]=0 with j < i
     *
     * @param s the search to follow
     * @param parts the parts of the pattern
     * @param allMatches vector to add occurrences to
     * @param exactMatchPositions a vector corresponding to the ranges for the
     * exact matches of the parts
     */
    void doRecSearch(const Search& s, std::vector<Substring>& parts,
                     std::vector<FMOcc<FMPosSFR>>& allMatches,
                     const std::vector<FMPosSFR>& exactMatchPositions,
                     const std::vector<std::vector<uint32_t>>& nodePaths) const;

    /**
     * @brief Match the pattern approximately in a naive way. All matches are at
     * most a certain edit distance away from the pattern. This function also
     * filters out the redundant matches using the maximal allowed edit
     * distance. Matching is done in a strain-free way.
     *
     * @param pattern the pattern to match
     * @param maxED the maximum edit distance
     * @return std::vector<FMOccSFR> - a vector with the resulting
     * occurrences in the compressed de Bruijn graph
     */
    std::vector<FMOccSFR> approxMatchesNaiveSFR(const std::string& pattern,
                                                length_t maxED) const;

    /**
     * @brief Splits the pattern into numParts parts, either by uniform range or
     * uniform size
     *
     * @param pattern the pattern to be split
     * @param parts the vector containing the substrings of this pattern,
     * will be cleared and filled during the execution of this method. If
     * the splitting fails for some reason, the vector will be empty
     * @param numparts, how many parts are needed
     * @param maxScore the maximum allowed edit distance
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partitionSFR(const std::string& pattern, std::vector<Substring>& parts,
                      const int& numParts, const int& maxScore,
                      std::vector<FMPosSFR>& exactMatchPositions,
                      std::vector<std::vector<uint32_t>>& nodePaths) const;

    /**
     * @brief Splits the pattern into numParts parts, such that each part has
     * the same size
     *
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numparts, how many parts are needed
     * @param maxScore, the maximum allowed edit distance
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void
    partitionUniformSFR(const std::string& pattern,
                        std::vector<Substring>& parts, const int& numParts,
                        const int& maxScore,
                        std::vector<FMPosSFR>& exactMatchPositions,
                        std::vector<std::vector<uint32_t>>& nodePaths) const;

    /**
     * @brief Splits the pattern into numParts parts, such that each search
     * carries the same weight (on average)
     *
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numparts, how many parts are needed
     * @param maxScore, the maximum allowed edit distance
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partitionOptimalStaticSFR(
        const std::string& pattern, std::vector<Substring>& parts,
        const int& numParts, const int& maxScore,
        std::vector<FMPosSFR>& exactMatchPositions,
        std::vector<std::vector<uint32_t>>& nodePaths) const;

    /**
     * @brief Splits the pattern into numParts parts, such that each part has
     * (approximately) the same range. The exactMatchRanges are also
     * calculated.
     *
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numparts, how many parts are needed
     * @param maxScore, the maximum allowed edit distance
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void
    partitionDynamicSFR(const std::string& pattern,
                        std::vector<Substring>& parts, const int& numParts,
                        const int& maxScore,
                        std::vector<FMPosSFR>& exactMatchPositions,
                        std::vector<std::vector<uint32_t>>& nodePaths) const;
};