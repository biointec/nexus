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

#include "strainfreemapper.h"

using namespace std;
// Check if this stays up to date with Luca's code
std::vector<FMOccSFR>
StrainFreeMapper::matchApproxSFR(const std::string& pattern,
                                 length_t maxED) const {

    strategy->index.resetCounters();

    if (maxED == 0) {

        SARangePair startRange = strategy->index.getCompleteRange();
        FMPosSFR pos(startRange, 0, 0);
        strategy->index.setDirection(FORWARD);
        Substring s(pattern, FORWARD);
        vector<uint32_t> nodePathRight;
        vector<uint32_t> nodePathLeft;
        pos.matchStringBidirectionally(s, nodePathLeft, nodePathRight);
        std::vector<FMOccSFR> result = {};
        if (pos.isValid()) {
            pos.setNodePath(nodePathRight);
            if (pos.getNodePath().empty()) {
                vector<uint32_t> nodePath;
                strategy->index.findNodeUnderK(pos, 0, nodePath);
                pos.setNodePath(nodePath);
            }
            FMOccSFR occ(pos, 0);
            result.emplace_back(occ);
        }
        return result;
    }
    // create the parts of the pattern
    vector<Substring> parts;

    // calculate how many parts there will be
    int numParts = strategy->calculateNumParts(maxED);
    // create the searches
    const vector<Search>& searches = strategy->createSearches(maxED);

    //  vector<SARangePair> exactMatchRanges(numParts);
    vector<FMPosSFR> exactMatchPositions(numParts, FMPosSFR());
    vector<vector<uint32_t>> nodePaths(numParts);
    partitionSFR(pattern, parts, numParts, maxED, exactMatchPositions,
                 nodePaths);

    if (parts.empty() || numParts * maxED >= pattern.size()) {
        if (!(strategy->getName() == "Naive backtracking")) {

            // splitting up was not viable just search the entire pattern
            cerr << "Warning: Normal bidirectional search was used as "
                    "entered pattern is too short "
                 << pattern.size() << endl;
        }

        return approxMatchesNaiveSFR(pattern, maxED);
    }

    // the  vector containing all matches in the suffix array
    vector<FMOcc<FMPosSFR>> allMatches;

    strategy->index.reserveStacks(numParts, pattern.length());
    // create the bit-parallel alignment matrices
    strategy->index.resetMatrices(
        parts.size()); // reset the alignment matrix that will
                       // be (possibly) used for each part

    // do all searches
    for (const Search& s : searches) {
        doRecSearch(s, parts, allMatches, exactMatchPositions, nodePaths);
    }

    // return all matches mapped to the text
    return strategy->index.filterStrainFreeMatches(allMatches, maxED);
}

void StrainFreeMapper::doRecSearch(
    const Search& s, vector<Substring>& parts,
    vector<FMOcc<FMPosSFR>>& allMatches,
    const vector<FMPosSFR>& exactMatchPositions,
    const std::vector<vector<uint32_t>>& nodePaths) const {

    if (s.getUpperBound(0) > 0) {
        // first part is allowed an error so start with an empty match
        s.setDirectionsInParts(parts);

        SARangePair startRange = strategy->index.getCompleteRange();
        FMOcc<FMPosSFR> startMatch = FMOcc<FMPosSFR>(startRange, 0, 0);
        vector<uint32_t> rightPath, leftPath;
        (strategy->*(strategy->startIdxPtrStrainFree))(
            s, startMatch, allMatches, parts, 0, leftPath, rightPath);
        return;
    }

    // first get the bidirectional match of first part
    int first = s.getPart(0);
    FMPosSFR startPos = exactMatchPositions[first];
    auto rightPath = nodePaths[first];

    vector<uint32_t> leftPath;

    if (startPos.isValid()) {
        s.setDirectionsInParts(parts);

        int partInSearch = 1;
        length_t exactLength = parts[first].size();

        while (s.getUpperBound(partInSearch) == 0) {
            // extend the exact match
            strategy->index.setDirection(s.getDirection(partInSearch - 1));

            startPos.matchStringBidirectionally(parts[s.getPart(partInSearch)],
                                                leftPath, rightPath);
            if (!startPos.isValid()) {
                return;
            }
            exactLength += parts[s.getPart(partInSearch)].size();
            partInSearch++;
        }

        startPos.setDepth(exactLength);
        startPos.setTrueDepth(exactLength);

        FMOcc<FMPosSFR> startMatch = FMOcc<FMPosSFR>(startPos, 0);

        (strategy->*(strategy->startIdxPtrStrainFree))(
            s, startMatch, allMatches, parts, partInSearch, leftPath,
            rightPath);
    }
}

std::vector<FMOccSFR>
StrainFreeMapper::approxMatchesNaiveSFR(const std::string& pattern,
                                        length_t maxED) const {

    vector<FMOcc<FMPosSFR>> occurrences =
        strategy->index.approxMatchesNaiveIntermediate(pattern, maxED);

    return strategy->index.filterStrainFreeMatches(occurrences, maxED);
}

void StrainFreeMapper::partitionSFR(const string& pattern,
                                    vector<Substring>& parts,
                                    const int& numParts, const int& maxScore,
                                    vector<FMPosSFR>& exactMatchPositions,
                                    vector<vector<uint32_t>>& nodePaths) const {

    parts.clear();

    if (numParts >= (int)pattern.size() || numParts == 1) {
        // no need of splitting up since all parts would be one
        // character or less or there is only one part
        return;
    }

    (this->*partitionPtrSFR)(pattern, parts, numParts, maxScore,
                             exactMatchPositions, nodePaths);
}

// Uniform Partitioning
void StrainFreeMapper::partitionUniformSFR(
    const string& pattern, vector<Substring>& parts, const int& numParts,
    const int& maxScore, vector<FMPosSFR>& exactMatchPositions,
    vector<vector<uint32_t>>& nodePaths) const {

    for (int i = 0; i < numParts; i++) {
        parts.emplace_back(pattern, (i * 1.0 / numParts) * pattern.size(),
                           ((i + 1) * 1.0 / numParts) * pattern.size());
    }
    // set end of final part correct
    parts.back().setEnd(pattern.size());

    // match the exactRanges for each part
    strategy->index.setDirection(FORWARD);

    for (int i = 0; i < numParts; i++) {
        vector<uint32_t> leftPath;
        exactMatchPositions[i].setRanges(strategy->index.getCompleteRange());
        exactMatchPositions[i].matchStringBidirectionally(parts[i], leftPath,
                                                          nodePaths[i]);
    }
}

// Static Partitioning
void StrainFreeMapper::partitionOptimalStaticSFR(
    const string& pattern, vector<Substring>& parts, const int& numParts,
    const int& maxScore, vector<FMPosSFR>& exactMatchPositions,
    vector<vector<uint32_t>>& nodePaths) const {

    strategy->setParts(pattern, parts, numParts, maxScore);

    // match the exactRanges for each part
    strategy->index.setDirection(FORWARD);

    for (int i = 0; i < numParts; i++) {
        vector<uint32_t> leftPath;
        exactMatchPositions[i].setRanges(strategy->index.getCompleteRange());
        exactMatchPositions[i].matchStringBidirectionally(parts[i], leftPath,
                                                          nodePaths[i]);
    }
}

// Dynamic Partitioning
void StrainFreeMapper::partitionDynamicSFR(
    const string& pattern, vector<Substring>& parts, const int& numParts,
    const int& maxScore, vector<FMPosSFR>& exactMatchPositions,
    vector<vector<uint32_t>>& nodePaths) const {
    throw runtime_error("Dynamic partitioning cannot be used with "
                        "strain-free matching.");
}