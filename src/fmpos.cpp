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

#include "fmpos.h"
#include "fmindexDBG.h"

using namespace std;

template <class positionClass>
thread_local ExtraCharPtr<positionClass> FMIndex<positionClass>::extraChar;
template <class positionClass>
thread_local length_t FMIndex<positionClass>::nodeCounter;

// ============================================================================
// CLASS FMPos
// ============================================================================

const void* FMPos::index;

void FMPos::extendFMPos(vector<FMPosExt<FMPos>>& stack) const {
    // Call the regular extendFMPos function from the index
    ((FMIndexDBG<FMPos>*)index)
        ->extendFMPos(getRanges(), stack, getDepth(), getTrueDepth());
}

const bool FMPos::separationIsNext() const {
    // Iterate over all separation characters
    for (int i = 0;
         i < ((FMIndexDBG<FMPos>*)index)->numberOfSeparationCharacters; i++) {
        SARangePair pairForNewChar;
        // Check if the current match can be extended using this separation
        // character
        if ((((FMIndexDBG<FMPos>*)index)->*FMIndex<FMPos>::extraChar)(
                i, getRanges(), pairForNewChar)) {
            return true;
        }
    }
    return false;
}

bool FMPos::compare(const FMOcc<FMPos>& rhs, int distance,
                    length_t shift) const {
    if (getRanges().getRangeSA().getBegin() !=
        rhs.getRanges().getRangeSA().getBegin()) {
        // the positions do not beging at the same position, first position is
        // smaller
        return getRanges().getRangeSA().getBegin() <
               rhs.getRanges().getRangeSA().getBegin();
    } else if (distance != rhs.getDistance()) {
        // begin is equal, better ed is smarter
        return distance < rhs.getDistance();
    } else if (getRanges().width() != rhs.getRanges().width()) {
        // shorter read is smaller...
        return getRanges().width() < rhs.getRanges().width();
    } else {
        // prefer no shift
        return shift < rhs.getShift();
    }
}

ostream& operator<<(ostream& o, const FMPos& pos) {
    return o << "SARange: " << pos.getRanges().getRangeSA()
             << "\tReverse SARange: " << pos.getRanges().getRangeSARev()
             << "\tDepth: " << pos.getDepth()
             << "\tTrue depth: " << pos.getTrueDepth();
}

// ============================================================================
// CLASS FMPosSFR
// ============================================================================

// Function pointers
ExtendFMPosAboveKPtr FMPosSFR::ExtendFMPosAboveK;
AddCharAboveKPtr FMPosSFR::AddCharAboveK;

void FMPosSFR::changeNodeRepresentation() {
    // Get the left boundary of the SA interval of the match
    length_t lb = ranges.getRangeSA().getBegin();
    // Find the first node of the path
    uint32_t id, l;
    ((FMIndexDBG<FMPosSFR>*)index)->findID(lb, id, l);
    // Add the node to the path
    appendToNodePath(id);
    // this is also the final node on the left (no left extension has happened
    // yet)
    finalNodeLeft = id;
    // Set the distances from the match to the ends of the node to their correct
    // starting values.
    setDistanceFromLeftEnd(l);
    setDistanceFromRightEnd(((FMIndexDBG<FMPosSFR>*)index)->G[id].len - l -
                            k_DBG);
}

void FMPosSFR::extendFMPos(vector<FMPosExt<FMPosSFR>>& stack) const {
    if (getTrueDepth() < k_DBG) {
        // The position does not have a node path yet
        // iterate over the entire alphabet, excluding the separation characters
        for (length_t i = 2; i < ((FMIndexDBG<FMPosSFR>*)index)->sigma.size();
             i++) {
            // Call the regular extending procedure on the index
            if (((FMIndexDBG<FMPosSFR>*)index)
                    ->extendFMPosIntermediary(getRanges(), stack, getDepth(), i,
                                              getTrueDepth()) &&
                getTrueDepth() == k_DBG - 1) {
                // If the position could be extended by the current character
                // and the new true depth equals the de Bruijn k parameter, the
                // node representation must be changed. In other words, the
                // corresponding node in the graph must be found.
                stack.back().changeNodeRepresentation();
            }
        }
        return;
    }

    // Call the extend procedure for positions that do have a node path. Which
    // function is called, depends on the current direction.
    (this->*ExtendFMPosAboveK)(stack);
}

void FMPosSFR::extendFMPosAboveKBackward(
    vector<FMPosExt<FMPosSFR>>& stack) const {

    // Get the first node of the current node path
    Node node = ((FMIndexDBG<FMPosSFR>*)index)->G[finalNodeLeft];
    // Get the index of the first suffix in the SA range of this position
    length_t indexInSA = getRanges().getRangeSA().getBegin();

    // Check whether we must switch to a new front node or not
    if (indexInSA != node.left_kmer_forward) {
        // We can add a character whilst staying in the current front node
        // Get the character that was added
        char c = ((FMIndexDBG<FMPosSFR>*)index)->bwt[indexInSA];
        // Check that it is no separation character
        if (!(c == '$' || c == '%')) {
            // Push this position onto the stack, along with the character that
            // was added.

            stack.emplace_back(c, *this);
            // Benchmarking
            ((FMIndexDBG<FMPosSFR>*)index)->nodeCounter++;
            // Get the added position
            FMPosExt<FMPosSFR>& newPosition = stack.back();
            // Increment the depth of the added position
            newPosition.incrementDepth();
            // Increment the true depth of the added position
            newPosition.incrementTrueDepth();
            // Get the index of the extended match in the SA by using the LF
            // property. Only one extension is possible, since we stay in the
            // same node.
            length_t newIndexInSA =
                ((FMIndexDBG<FMPosSFR>*)index)->findLF(indexInSA, false);
            // Update the ranges of the added position. Only the range in the SA
            // changes.
            newPosition.setRanges(SARangePair(
                Range(newIndexInSA, newIndexInSA + node.multiplicity),
                getRanges().getRangeSARev()));
            // Update the distance of the start of the match to the start of the
            // first node
            newPosition.setDistanceFromLeftEnd(distanceFromLeftEnd - 1);
        }
    } else {
        // We cannot add a character whilst staying in the current front node.
        // Hence, we need to investigate all predecessors of the front node.

        // Get the identifier of the first node of the path
        uint32_t id = finalNodeLeft;
        // Create a variable in which the id of the predecessor will be stored
        uint32_t id_predecessor;
        // Iterate over all possible characters for extension. Separation
        // characters are excluded.
        for (size_t i = 2; i < ((FMIndexDBG<FMPosSFR>*)index)->sigma.size();
             i++) {
            // Check if there exists a predecessor that is the result of
            // extension with the current character
            if (((FMIndexDBG<FMPosSFR>*)index)
                    ->jumpToPredecessorWithChar(id, id_predecessor, i)) {
                // A predecessor exists
                // Push this position onto the stack, along with the character
                // that was added.
                stack.emplace_back(((FMIndexDBG<FMPosSFR>*)index)->sigma.i2c(i),
                                   *this);
                // Benchmarking
                ((FMIndexDBG<FMPosSFR>*)index)->nodeCounter++;
                // Get the added position
                FMPosExt<FMPosSFR>& newPosition = stack.back();
                // Increment the depth of the added position
                newPosition.incrementDepth();
                // Increment the true depth of the added position
                newPosition.incrementTrueDepth();
                // Prepend the predecessor to the node path
                newPosition.prependToNodePath(id_predecessor);
                // Get the node attributes corresponding to the predecessor
                Node newNode =
                    ((FMIndexDBG<FMPosSFR>*)index)->G[id_predecessor];
                // Update the ranges of the added position. Only the range in
                // the SA changes. It can be deduced from the attributes of the
                // predecessor.
                newPosition.setRanges(SARangePair(
                    Range(newNode.right_kmer_forward,
                          newNode.right_kmer_forward + newNode.multiplicity),
                    getRanges().getRangeSARev()));
                // Update the distance of the start of the match to the start of
                // the first node
                newPosition.setDistanceFromLeftEnd(
                    ((FMIndexDBG<FMPosSFR>*)index)->G[id_predecessor].len -
                    k_DBG);
            }
        }
    }
}

void FMPosSFR::extendFMPosAboveKForward(
    vector<FMPosExt<FMPosSFR>>& stack) const {
    // Get the last node of the current node path
    Node node = ((FMIndexDBG<FMPosSFR>*)index)->G[finalNodeRight];
    // Get the index of the first suffix in the reverse SA range of this
    // position
    length_t indexInRevSA = getRanges().getRangeSARev().getBegin();
    // Check whether we must switch to a new end node or not
    if (indexInRevSA != node.right_kmer_reverse) {
        // We can add a character whilst staying in the current end node
        // Get the character that was added
        char c = ((FMIndexDBG<FMPosSFR>*)index)->revbwt[indexInRevSA];
        // Check that it is no separation character
        if (!(c == '$' || c == '%')) {
            // Push this position onto the stack, along with the character that
            // was added.
            stack.emplace_back(c, *this);
            // Benchmarking
            ((FMIndexDBG<FMPosSFR>*)index)->nodeCounter++;
            // Get the added position
            FMPosExt<FMPosSFR>& newPosition = stack.back();
            // Increment the depth of the added position
            newPosition.incrementDepth();
            // Increment the true depth of the added position
            newPosition.incrementTrueDepth();
            // Get the index of the extended match in the reverse SA by using
            // the LF property on the reverse BWT. Only one extension is
            // possible, since we stay in the same node.
            length_t newIndexInRevSA =
                ((FMIndexDBG<FMPosSFR>*)index)->findLF(indexInRevSA, true);
            // Update the ranges of the added position. Only the range in the
            // reverse SA changes.
            newPosition.setRanges(SARangePair(
                getRanges().getRangeSA(),
                Range(newIndexInRevSA, newIndexInRevSA + node.multiplicity)));
            // Update the distance of the end of the match to the end of the
            // last node
            newPosition.setDistanceFromRightEnd(distanceFromRightEnd - 1);
        }
    } else {
        // We cannot add a character whilst staying in the current end node.
        // Hence, we need to investigate all successors of the end node.

        // Get the identifier of the last node of the path
        uint32_t id = finalNodeRight;
        // Create a variable in which the id of the successor will be stored
        uint32_t id_successor;
        // Iterate over all possible characters for extension. Separation
        // characters are excluded.
        for (size_t i = 2; i < ((FMIndexDBG<FMPosSFR>*)index)->sigma.size();
             i++) {
            // Check if there exists a successor that is the result of
            // extension with the current character
            if (((FMIndexDBG<FMPosSFR>*)index)
                    ->jumpToSuccessorWithChar(id, id_successor, i)) {
                // A successor exists
                // Push this position onto the stack, along with the character
                // that was added.
                stack.emplace_back(((FMIndexDBG<FMPosSFR>*)index)->sigma.i2c(i),
                                   *this);
                // Benchmarking
                ((FMIndexDBG<FMPosSFR>*)index)->nodeCounter++;
                // Get the added position
                FMPosExt<FMPosSFR>& newPosition = stack.back();
                // Increment the depth of the added position
                newPosition.incrementDepth();
                // Increment the true depth of the added position
                newPosition.incrementTrueDepth();
                // Prepend the successor to the node path
                newPosition.appendToNodePath(id_successor);
                // Get the node attributes corresponding to the successor
                Node newNode = ((FMIndexDBG<FMPosSFR>*)index)->G[id_successor];
                // Update the ranges of the added position. Only the range in
                // the reverse SA changes. It can be deduced from the attributes
                // of the successor.
                newPosition.setRanges(SARangePair(
                    getRanges().getRangeSA(),
                    Range(newNode.left_kmer_reverse,
                          newNode.left_kmer_reverse + newNode.multiplicity)));
                // Update the distance of the end of the match to the end of the
                // last node
                newPosition.setDistanceFromRightEnd(
                    ((FMIndexDBG<FMPosSFR>*)index)->G[id_successor].len -
                    k_DBG);
            }
        }
    }
}

bool FMPosSFR::addCharAboveKBackward(const char& c) {
    // Get the first node of the current node path
    Node node = ((FMIndexDBG<FMPosSFR>*)index)->G[finalNodeLeft];
    // Get the index of the first suffix in the SA range of this position
    length_t indexInSA = getRanges().getRangeSA().getBegin();
    // Check whether we must switch to a new front node or not
    if (indexInSA != node.left_kmer_forward) {
        // We can add a character whilst staying in the current front node
        // Check if the next character in the current front node is equal to the
        // character c we want to add
        if (((FMIndexDBG<FMPosSFR>*)index)->bwt[indexInSA] != c) {
            // The characters are not equal, so character c cannot be added. To
            // indicate that this position is invalid, we set the ranges to
            // empty.
            setRanges(((FMIndexDBG<FMPosSFR>*)index)->getEmptyRange());
            return false;
        }
        // Increment the depth
        incrementDepth();
        // Increment the true depth
        incrementTrueDepth();
        // Get the index of the extended match in the SA by using the LF
        // property. Only one extension is possible, since we stay in the same
        // node.
        length_t newIndexInSA =
            ((FMIndexDBG<FMPosSFR>*)index)->findLF(indexInSA, false);
        // Update the ranges. Only the range in the SA changes.
        setRanges(
            SARangePair(Range(newIndexInSA, newIndexInSA + node.multiplicity),
                        getRanges().getRangeSARev()));
        // Update the distance of the start of the match to the start of the
        // first node
        setDistanceFromLeftEnd(distanceFromLeftEnd - 1);
        return true;
    } else {
        // We cannot add a character whilst staying in the current front node.
        // Hence, we need to investigate the possible predecessors of the front
        // node.

        // Get the index in alfabet sigma that corresponds to character c
        int posInAlphabet =
            ((FMIndexDBG<FMPosSFR>*)index)->sigma.c2i((unsigned char)c);
        // Check that this index is valid
        if (posInAlphabet > -1) {
            // Get the identifier of the first node of the path
            uint32_t id = finalNodeLeft;
            // Create a variable in which the id of the predecessor will be
            // stored
            uint32_t id_predecessor;
            // Check if there exists a predecessor that is the result of
            // extension with character c
            if (((FMIndexDBG<FMPosSFR>*)index)
                    ->jumpToPredecessorWithChar(id, id_predecessor,
                                                posInAlphabet)) {
                // A predecessor exists
                // Increment the depth
                incrementDepth();
                // Increment the true depth
                incrementTrueDepth();
                // Prepend the predecessor to the node path
                prependToNodePath(id_predecessor);
                // Get the node attributes corresponding to the predecessor
                Node newNode =
                    ((FMIndexDBG<FMPosSFR>*)index)->G[id_predecessor];
                // Update the ranges. Only the range in the SA changes. It can
                // be deduced from the attributes of the predecessor.
                setRanges(SARangePair(
                    Range(newNode.right_kmer_forward,
                          newNode.right_kmer_forward + newNode.multiplicity),
                    getRanges().getRangeSARev()));
                // Update the distance of the start of the match to the start of
                // the first node
                setDistanceFromLeftEnd(
                    ((FMIndexDBG<FMPosSFR>*)index)->G[id_predecessor].len -
                    k_DBG);
                return true;
            }
        }
        // The index in alfabet sigma is not valid or there exists no
        // predecessor that is the result of the extension of the current match
        // with character c. To indicate that this position is invalid, we set
        // the ranges to empty.
        setRanges(((FMIndexDBG<FMPosSFR>*)index)->getEmptyRange());
        return false;
    }
}

bool FMPosSFR::addCharAboveKForward(const char& c) {
    // Get the last node of the current node path
    Node node = ((FMIndexDBG<FMPosSFR>*)index)->G[finalNodeRight];
    // Get the index of the first suffix in the reverse SA range of this
    // position
    length_t indexInRevSA = getRanges().getRangeSARev().getBegin();
    // Check whether we must switch to a new end node or not
    if (indexInRevSA != node.right_kmer_reverse) {
        // We can add a character whilst staying in the current end node
        // Check if the next character in the current end node is equal to the
        // character c we want to add
        if (((FMIndexDBG<FMPosSFR>*)index)->revbwt[indexInRevSA] != c) {
            // The characters are not equal, so character c cannot be added. To
            // indicate that this position is invalid, we set the ranges to
            // empty.
            setRanges(((FMIndexDBG<FMPosSFR>*)index)->getEmptyRange());
            return false;
        }
        // Increment the depth
        incrementDepth();
        // Increment the true depth
        incrementTrueDepth();
        // Get the index of the extended match in the reverse SA by using the LF
        // property on the reverse BWT. Only one extension is possible, since we
        // stay in the same node.
        length_t newIndexInRevSA =
            ((FMIndexDBG<FMPosSFR>*)index)->findLF(indexInRevSA, true);
        // Update the ranges. Only the range in the reverse SA changes.
        setRanges(SARangePair(
            getRanges().getRangeSA(),
            Range(newIndexInRevSA, newIndexInRevSA + node.multiplicity)));
        // Update the distance of the end of the match to the end of the
        // last node
        setDistanceFromRightEnd(distanceFromRightEnd - 1);
        return true;
    } else {
        // We cannot add a character whilst staying in the current end node.
        // Hence, we need to investigate the possible successors of the end
        // node.

        // Get the index in alfabet sigma that corresponds to character c
        int posInAlphabet =
            ((FMIndexDBG<FMPosSFR>*)index)->sigma.c2i((unsigned char)c);
        // Check that this index is valid
        if (posInAlphabet > -1) {
            // Get the identifier of the last node of the path
            uint32_t id = finalNodeRight;
            // Create a variable in which the id of the successor will be
            // stored
            uint32_t id_successor;
            // Check if there exists a successor that is the result of
            // extension with character c
            if (((FMIndexDBG<FMPosSFR>*)index)
                    ->jumpToSuccessorWithChar(id, id_successor,
                                              posInAlphabet)) {
                // A successor exists
                // Increment the depth
                incrementDepth();
                // Increment the true depth
                incrementTrueDepth();
                // Append the successor to the node path
                appendToNodePath(id_successor);
                // Get the node attributes corresponding to the successor
                Node newNode = ((FMIndexDBG<FMPosSFR>*)index)->G[id_successor];
                // Update the ranges. Only the range in the reverse SA changes.
                // It can be deduced from the attributes of the successor.
                setRanges(SARangePair(
                    getRanges().getRangeSA(),
                    Range(newNode.left_kmer_reverse,
                          newNode.left_kmer_reverse + newNode.multiplicity)));
                // Update the distance of the end of the match to the end of the
                // last node
                setDistanceFromRightEnd(
                    ((FMIndexDBG<FMPosSFR>*)index)->G[id_successor].len -
                    k_DBG);
                return true;
            }
        }
        // The index in alfabet sigma is not valid or there exists no
        // successor that is the result of the extension of the current match
        // with character c. To indicate that this position is invalid, we set
        // the ranges to empty.
        setRanges(((FMIndexDBG<FMPosSFR>*)index)->getEmptyRange());
        return false;
    }
}

const bool FMPosSFR::separationIsNext() const {
    // First check which node representation we are in. This can be done by
    // checking whether the node path is empty.
    if (!hasNodePath()) {
        // The position does not have a node path yet

        // Iterate over all separation characters
        for (int i = 0;
             i < ((FMIndexDBG<FMPosSFR>*)index)->numberOfSeparationCharacters;
             i++) {
            SARangePair pairForNewChar;
            // Check if the current match can be extended using this
            // separation character. This is done by trying to add a
            // character based on the underlying bidirectional FM-index.
            if ((((FMIndexDBG<FMPosSFR>*)index)
                     ->*FMIndexDBG<FMPosSFR>::extraChar)(i, getRanges(),
                                                         pairForNewChar)) {
                return true;
            }
        }
        return false;
    }
    // The position does have a node path

    // Iterate over all separation characters
    for (int i = 0;
         i < ((FMIndexDBG<FMPosSFR>*)index)->numberOfSeparationCharacters;
         i++) {
        // Find the separation character that corresponds to index i in alfabet
        // sigma.
        unsigned char c = ((FMIndexDBG<FMPosSFR>*)index)->sigma.i2c(i);
        // Create a copy of the current position to execute the character adding
        // operation on (this object will be altered)
        FMPosSFR copy = *this;
        // Check if the current match can be extended using this separation
        // character. This is done by trying to add a character based on the
        // compressed de Bruijn graph.
        if ((copy.*AddCharAboveK)(c)) {
            return true;
        }
    }
    return false;
}

bool FMPosSFR::compare(const FMOcc<FMPosSFR>& rhs, int distance,
                       length_t shift) const {

    const auto& pos = rhs.getPosition();
    const auto& lhsNodePath = getNodePath();
    const auto& rhsNodePath = pos.getNodePath();
    if (lhsNodePath != rhsNodePath) {
        // The node paths are not identical, they are compared elementwise to
        // see which is smaller
        return lhsNodePath < rhsNodePath;
    } else if (distanceFromLeftEnd !=
               rhs.getPosition().getDistanceFromLeftEnd()) {
        // The distances from the start of the match to the start of the first
        // node are not identical, smaller distance wins
        return distanceFromLeftEnd < rhs.getPosition().getDistanceFromLeftEnd();
    } else if (distance != rhs.getDistance()) {
        // begin is equal, better ed is smarter
        return distance < rhs.getDistance();
    } else if (trueDepth != pos.getTrueDepth()) {
        // begin is equal, shorter match is smaller
        return trueDepth < pos.getTrueDepth();
    } else {
        // prefer no shift
        return shift < rhs.getShift();
    }
}

void FMPosSFR::matchStringBidirectionally(Substring& pattern,
                                          vector<uint32_t>& nodePathLeft,
                                          vector<uint32_t>& nodePathRight) {

    // Get the range pair corresponding to this position
    SARangePair rangePair = getRanges();

    // Calculate the number of positions that still need to be matched until a
    // node path can be found.
    size_t bound = (k_DBG > trueDepth) ? k_DBG - trueDepth : 0;

    // Match characters using the underlying bidirectional FM-index until either
    // the end of the pattern is reached, or the node representation must
    // change.
    for (length_t i = 0; i < min(pattern.size(), (length_t)bound); i++) {
        // Get the next character from the pattern
        char c = pattern[i];
        // Check if the current match can be extended using this next character.
        // This is done by trying to add the character based on the underlying
        // bidirectional FM-index.
        if (!((FMIndexDBG<FMPosSFR>*)index)->addChar(c, rangePair)) {
            // The character cannot be matched. Update the ranges such that they
            // are empty.
            setRanges(rangePair);
            return;
        }
        // Increment the depth
        incrementDepth();
        // Increment the true depth
        incrementTrueDepth();
    }

    // Update the range pair to the first part of the exact match
    setRanges(rangePair);

    if (pattern.size() < bound) {
        // We have reached the end of the pattern and thus do not need to
        // continue matching based on the compressed de Bruijn graph.
        return;
    }

    if (trueDepth == k_DBG && !hasNodePath()) {
        // The total true depth of the current match equals the de Bruijn
        // parameter k and the node representation must thus be changed. In
        // other words, the node path is found.
        changeNodeRepresentation();
        nodePathRight.emplace_back(finalNodeRight);
    }

    // Match characters using the compressed de Bruijn graph until the end of
    // the pattern is reached.
    for (length_t i = bound; i < pattern.size(); i++) {

        // Get the next character from the pattern
        char c = pattern[i];

        // Check if the current match can be extended using this character. This
        // is done by trying to add the character based on the compressed de
        // Bruijn graph.
        if (!(this->*AddCharAboveK)(c)) {
            return;
        }

        // extend the node paths (right always has a nodepath, left maybe)
        nodePathRight.resize(numberOfNodesRight);
        nodePathRight.back() = finalNodeRight;
        nodePathLeft.resize(numberOfNodesLeft);
        if (nodePathLeft.size()) {
            nodePathLeft.back() = finalNodeLeft;
        }
        ((FMIndexDBG<FMPosSFR>*)index)->nodeCounter++;
    }
}

ostream& operator<<(ostream& o, const FMPosSFR& pos) {
    o << "SARange: " << pos.getRanges().getRangeSA()
      << "\tReverse SARange: " << pos.getRanges().getRangeSARev()
      << "\tDepth: " << pos.getDepth() << "\tTrue depth: " << pos.getTrueDepth()
      << "\tLeftmost node: " << pos.getFinalNodeLeft()
      << "\tnumber of nodes left" << pos.getNumberOfNodesLeft()
      << "\tRightmost node: " << pos.getFinalNodeRight()
      << "\tnumber of nodes left" << pos.getNumberOfNodesRight()
      << "\tNode path: " << pos.getNodePath()[0];
    for (length_t i = 1; i < pos.getNodePath().size(); i++) {
        o << "," << pos.getNodePath()[i];
    }

    return o;
}