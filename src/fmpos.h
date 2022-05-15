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

#include "sarangepair.h"
#include "substring.h"

#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>

// ============================================================================
// CLASS FMPos
// ============================================================================

template <class positionClass> class FMIndex;
template <class positionClass> class FMIndexDBG;
template <class positionClass> class FMPosExt;
template <class positionClass> class FMOcc;
class FMPosSFR;

/**
 * A position in the bidirectional FM-index.
 */
class FMPos {
  protected:
    // the ranges over the suffix arrays
    SARangePair ranges;
    // the depth of the prefix of the suffixes of this position with respect to
    // the current score matrix
    length_t depth;
    // the true depth of the total current match
    length_t trueDepth;
    // a pointer to the index that is being used in this session
    static const void* index;

  public:
    /**
     * @brief Construct a new FMPos object. Default constructor for empty
     * position (= empty ranges and depth of zero).
     *
     */
    FMPos() : ranges(SARangePair()), depth(0), trueDepth(0) {
    }

    /**
     * @brief Construct a new FMPos object.
     *
     * @param ranges the ranges over the suffix arrays
     * @param depth the depth of the prefix of the suffixes of this position
     * with respect to the current score matrix
     * @param trueDepth the true depth of the total current match. Default = -1.
     */
    FMPos(SARangePair& ranges, int depth, int trueDepth = -1)
        : ranges(ranges), depth(depth), trueDepth(trueDepth) {
        if (trueDepth == -1) {
            this->trueDepth = depth;
        }
    }

    /**
     * @brief Static function. Set the index for this session.
     *
     * @param _index - a pointer to the index that is being used in this session
     */
    static void setIndex(const void* _index) {
        index = _index;
    }

    /**
     * @brief Get the ranges over the suffix arrays
     *
     * @return const SARangePair& - the ranges over the suffix arrays
     */
    const SARangePair& getRanges() const {
        return ranges;
    }

    /**
     * @brief Get the depth of the prefix of the suffixes of this position with
     * respect to the current score matrix
     *
     * @return const length_t& - the depth of the prefix of the suffixes of this
     * position with respect to the current score matrix
     */
    const length_t& getDepth() const {
        return depth;
    }

    /**
     * @brief Get the true depth of the total current match
     *
     * @return const length_t - the true depth of the total current match
     */
    const length_t getTrueDepth() const {
        return trueDepth;
    }

    /**
     * @brief Set the ranges over the suffix arrays
     *
     * @param ranges the ranges over the suffix arrays
     */
    void setRanges(SARangePair ranges) {
        this->ranges = ranges;
    }

    /**
     * @brief Set the depth of the prefix of the suffixes of this position with
     * respect to the current score matrix
     *
     * @param depth - the depth of the prefix of the suffixes of this position
     * with respect to the current score matrix
     */
    void setDepth(length_t depth) {
        this->depth = depth;
    }

    /**
     * @brief Set the true depth of the total current match
     *
     * @param trueDepth the true depth of the total current match
     */
    void setTrueDepth(length_t trueDepth) {
        this->trueDepth = trueDepth;
    }

    /**
     * @brief Increment the depth of the prefix of the suffixes of this position
     * with respect to the current score matrix
     *
     */
    void incrementDepth() {
        depth++;
    }

    /**
     * @brief Increment the true depth of the total current match
     *
     */
    void incrementTrueDepth() {
        trueDepth++;
    }

    /**
     * @brief Check if the FMPos is valid. In other words, check if the ranges
     * are empty.
     *
     * @return true if the ranges are not empty
     * @return false otherwise
     */
    bool isValid() const {
        return !ranges.empty();
    }

    /**
     * @brief Get the size of the FMPos. In other words, get the number of times
     * that the text corresponding to this position occurs in the original text.
     *
     * @return length_t - the size of the FMPos
     */
    length_t size() const {
        return ranges.width();
    }

    virtual void reverseNodePath() {
        throw std::runtime_error(
            "The FMPos object has no node path to reverse.");
    }

    virtual void
    updateNodeStackWithNodePath(std::vector<uint32_t>& leftNodes,
                                std::vector<uint32_t>& rightNodes) const {
        throw std::runtime_error(
            "The FMPos object has no node path to update the stacks with");
    }

    /**
     * @brief Push all the children corresponding to the node with ranges
     * equal to ranges.
     *
     * @param stack the stack to push the children on
     */
    void extendFMPos(std::vector<FMPosExt<FMPos>>& stack) const;

    /**
     * @brief Check if one of the children of this position is a separation
     * character.
     *
     * @return true if one of the children of this position is a separation
     * character
     * @return false otherwise
     */
    virtual const bool separationIsNext() const;

    /**
     * @brief Operator overloading, two FMPos are equal if their ranges and
     * depth are equal
     *
     * @param rhs the FMPos to compare to this
     * @return true if equal
     * @return false if not equal
     */
    bool operator==(const FMPos& rhs) const {
        return ranges == rhs.getRanges() && depth == rhs.getDepth() &&
               trueDepth == rhs.getTrueDepth();
    }

    /**
     * @brief Compare two FMOcc objects. This function is defined here, in the
     * FMPos class, since a distinction needs to be made between FMOcc objects
     * containing a FMPos versus FMOcc object containing a FMPosSFR. This
     * function is called on the FMOcc object that contains this FMPos object.
     *
     * @param rhs The other FMOcc object to compare this one to
     * @param distance The distance corresponding to this FMOcc object
     * @param shift The shift corresponding to this FMOcc object
     * @return true if this FMOcc object is smaller than rhs
     * @return false otherwise.
     */
    bool compare(const FMOcc<FMPos>& rhs, int distance, length_t shift) const;

    /**
     * @brief Operator overloading. Outputs the FMPos object to the outputstream
     *
     * @param os the output stream
     * @param pos the FMPos to print
     * @return std::ostream& - the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const FMPos& pos);
};

// ============================================================================
// CLASS FMPosSFR
// ============================================================================

/**
 * A position in the compressed de Bruijn graph. This class is only used during
 * strain-free matching.
 */

using ExtendFMPosAboveKPtr =
    void (FMPosSFR::*)(std::vector<FMPosExt<FMPosSFR>>&) const;

using AddCharAboveKPtr = bool (FMPosSFR::*)(const char& c);

class FMPosSFR : public FMPos {
  protected:
    // Pointer to The node path in the graph that represents the current match
    std::shared_ptr<std::vector<uint32_t>> ptrToNodePath;

    uint32_t finalNodeLeft;  // the final node of the path on the left side
    uint32_t finalNodeRight; // the final node of the path on the right side
    uint32_t numberOfNodesLeft = 0;  // the number of nodes on the left side
    uint32_t numberOfNodesRight = 0; // the number of nodes on the right side

    // The distance from the start of the match to the start of the first node
    // of the node path
    uint32_t distanceFromLeftEnd;
    // The distance from the end of the match to the end of the last node of the
    // node path
    uint32_t distanceFromRightEnd;

    // Note that the SARangePair of the superclass does not necessarily
    // represent a contiguous match in the original text anymore. The SA range
    // represents the left bound of the match in the first node of the node
    // path. The reverse SA range represents the right bound of the match in the
    // last node of the node path.

    // pointer to extendFMPos method (for direction)
    static ExtendFMPosAboveKPtr ExtendFMPosAboveK;

    // pointer to extra char method (for direction)
    static AddCharAboveKPtr AddCharAboveK;

  public:
    /**
     * @brief Construct a new FMPosSFR object. Default constructor for empty
     * position (= empty ranges, depth, true depth and distances of zero).
     *
     */
    FMPosSFR()
        : FMPos(), ptrToNodePath(nullptr), numberOfNodesLeft(0),
          numberOfNodesRight(0), distanceFromLeftEnd(0),
          distanceFromRightEnd(0) {
    }

    /**
     * @brief Construct a new FMPosSFR object
     *
     * @param ranges the SA range represents the left bound of the match in the
     * first node of the node path. The reverse SA range represents the right
     * bound of the match in the last node of the node path.
     * @param depth the depth of the prefix of the suffixes of this position
     * with respect to the current score matrix
     * @param trueDepth the true depth of the total current match. Default = 0.
     * @param distanceFromLeftEnd The distance from the start of the match to
     * the start of the first node of the node path. Default = 0.
     * @param distanceFromRightEnd The distance from the end of the match to the
     * end of the last node of the node path. Default = 0.
     */
    FMPosSFR(SARangePair& ranges, int depth, int trueDepth = 0,
             int distanceFromLeftEnd = 0, int distanceFromRightEnd = 0)
        : FMPos(ranges, depth, trueDepth), ptrToNodePath(nullptr),
          numberOfNodesLeft(0), numberOfNodesRight(0),
          distanceFromLeftEnd(distanceFromLeftEnd),
          distanceFromRightEnd(distanceFromRightEnd) {
    }

    /**
     * @brief Get the distance from the start of the match to the start of the
     * first node of the node path
     *
     * @return const uint32_t& - the distance from the start of the match to the
     * start of the first node of the node path
     */
    const uint32_t& getDistanceFromLeftEnd() const {
        return distanceFromLeftEnd;
    }

    /**
     * @brief Set the distance from the start of the match to the start of the
     * first node of the node path
     *
     * @param distanceFromLeftEnd the new distance from the start of the match
     * to the start of the first node of the node path
     */
    void setDistanceFromLeftEnd(uint32_t distanceFromLeftEnd) {
        this->distanceFromLeftEnd = distanceFromLeftEnd;
    }

    /**
     * @brief Get the distance from the end of the match to the end of the last
     * node of the node path.
     *
     * @return const uint32_t& - the distance from the end of the match to the
     * end of the last node of the node path.
     */
    const uint32_t& getDistanceFromRightEnd() const {
        return distanceFromRightEnd;
    }

    /**
     * @brief Set the distance from the end of the match to the end of the last
     * node of the node path.
     *
     * @param distanceFromRightEnd the new distance from the end of the match to
     * the end of the last node of the node path.
     */
    void setDistanceFromRightEnd(uint32_t distanceFromRightEnd) {
        this->distanceFromRightEnd = distanceFromRightEnd;
    }

    /**
     * @brief Swap the distance from the start of the match to the start of the
     * first node of the node path and the distance from the end of the match to
     * the end of the last node of the node path.
     *
     */
    void swapDistances() {
        std::swap(distanceFromLeftEnd, distanceFromRightEnd);
    }

    /**
     * @brief Get the node path in the graph that represents the current match
     *
     * @return const std::vector<uint32_t>& - the node path in the graph that
     * represents the current match
     */
    const std::vector<uint32_t>& getNodePath() const {
        assert(ptrToNodePath.get() != nullptr);
        return *ptrToNodePath.get();
    }

    /**
     * @brief Set the node path in the graph that represents the current match
     *
     * @param nodePath the node path in the graph that represents the current
     * match
     */
    void setNodePath(std::vector<uint32_t>& nodePath) {
        ptrToNodePath = std::make_shared<std::vector<uint32_t>>(nodePath);
    }

    uint32_t getFinalNodeLeft() const {
        return finalNodeLeft;
    }

    uint32_t getFinalNodeRight() const {
        return finalNodeRight;
    }

    uint32_t getNumberOfNodesLeft() const {
        return numberOfNodesLeft;
    }

    uint32_t getNumberOfNodesRight() const {
        return numberOfNodesRight;
    }

    bool hasNodePath() const {
        return numberOfNodesLeft + numberOfNodesRight > 0;
    }
    virtual void reverseNodePath() override {
        std::reverse(ptrToNodePath.get()->begin(), ptrToNodePath.get()->end());
    }

    /**
     * Update the left and right nodes to correspond to the node path of this
     * position
     * @param leftNodes the path to the left
     * @param rightNodes the path to the right
     */
    virtual void updateNodeStackWithNodePath(
        std::vector<uint32_t>& leftNodes,
        std::vector<uint32_t>& rightNodes) const override {
        // reset the left and right nodes
        leftNodes.resize(getNumberOfNodesLeft());
        if (leftNodes.size() > 0) {
            leftNodes.back() = getFinalNodeLeft();
        }

        rightNodes.resize(getNumberOfNodesRight());
        if (rightNodes.size() > 0) {
            rightNodes.back() = getFinalNodeRight();
        }
    }

    /**
     * @brief Append a node to the end of the node path in the graph that
     * represents the current match
     *
     * @param nodeID the identifier of the node that must be appended
     */
    void appendToNodePath(int nodeID) {
        // nodePath.emplace_back(nodeID);
        //  append is to the right
        numberOfNodesRight++;
        finalNodeRight = nodeID;
    }

    /**
     * @brief Prepend a node to the front of the node path in the graph that
     * represents the current match
     *
     * @param nodeID the identifier of the node that must be prepended
     */
    void prependToNodePath(int nodeID) {
        // nodePath.insert(nodePath.begin(), nodeID);
        //  prepend is to the left
        numberOfNodesLeft++;
        finalNodeLeft = nodeID;
    }

    /**
     * @brief Operator overloading, two FMPosSFR objects are equal if their
     * ranges, depth, true depth and node path are equal.
     *
     * @param rhs the FMPosSFR object to compare to this
     * @return true if this is equal to rhs
     * @return false otherwise
     */
    bool operator==(const FMPosSFR& rhs) const {
        return ranges == rhs.getRanges() && depth == rhs.getDepth() &&
               trueDepth == rhs.getTrueDepth() &&
               getNodePath() == rhs.getNodePath();
    }

    /**
     * @brief Check if this FMPosSFR object is valid. In other words, check
     * that its ranges are not empty.
     *
     * @return true if the ranges are not empty
     * @return false otherwise
     */
    bool isValid() const {
        return !ranges.empty();
    }

    /**
     * @brief Change the node representation of this FMPosSFR object. This
     * function must be called when the true depth reaches the de Bruijn
     * parameter k. At this point, the node path in the graph can be found.
     *
     */
    void changeNodeRepresentation();

    /**
     * @brief Push all the children corresponding to the node with this node
     * path and ranges. Separation characters cannot be added.
     *
     * @param stack the stack to push the children on
     */
    void extendFMPos(std::vector<FMPosExt<FMPosSFR>>& stack) const;

    /**
     * @brief Push all the children corresponding to the node with this node
     * path and ranges. Specifically, we move backward or to the left in the
     * graph. This function is only called when the FMPosSFR object already
     * has a node path in the graph. Separation characters cannot be added.
     *
     * @param stack the stack to push the children on
     */
    void
    extendFMPosAboveKBackward(std::vector<FMPosExt<FMPosSFR>>& stack) const;

    /**
     * @brief Push all the children corresponding to the node with this node
     * path and ranges. Specifically, we move forward or to the right in the
     * graph. This function is only called when the FMPosSFR object already
     * has a node path in the graph. Separation characters cannot be added.
     *
     * @param stack the stack to push the children on
     */
    void extendFMPosAboveKForward(std::vector<FMPosExt<FMPosSFR>>& stack) const;

    /**
     * @brief Add one character to the front or left of the match and update
     * the FMPosSFR attributes. If the character cannot be added, the range
     * will be set to an empty range. This function is only called when the
     * FMPosSFR object already has a node path in the graph. Separation
     * characters can be added.
     *
     * @param c the character to be added to the front of the current match
     * @return true if the character can be added
     * @return false if the character cannot be added
     */
    bool addCharAboveKBackward(const char& c);

    /**
     * @brief Add one character to the back or right of the match and update
     * the FMPosSFR attributes. If the character cannot be added, the range
     * will be set to an empty range. This function is only called when the
     * FMPosSFR object already has a node path in the graph. Separation
     * characters can be added.
     *
     * @param c the character to be added to the back of the current match
     * @return true if the character can be added
     * @return false if the character cannot be added
     */
    bool addCharAboveKForward(const char& c);

    /**
     * @brief Check if one of the children of this position is a separation
     * character. Children are found in a strain-free manner.
     *
     * @return true if one of the children of this position is a separation
     * character
     * @return false otherwise
     */
    virtual const bool separationIsNext() const override;

    /**
     * @brief Compare two FMOcc objects. This function is defined here, in
     * the FMPosSFR class, since a distinction needs to be made between FMOcc
     * objects containing a FMPos versus FMOcc object containing a FMPosSFR.
     * This function is called on the FMOcc object that contains this
     * FMPosSFR object.
     *
     * @param rhs The other FMOcc object to compare this one to
     * @param distance The distance corresponding to this FMOcc object
     * @param shift The shift corresponding to this FMOcc object
     * @return true true if this FMOcc object is smaller than rhs
     * @return false false otherwise.
     */
    bool compare(const FMOcc<FMPosSFR>& rhs, int distance,
                 length_t shift) const;

    /**
     * @brief Match a pattern to the implicit de Bruijn graph in an exact
     * way. Whether or not a node path is found, depends on the true depth
     * of the match.
     *
     * @param pattern The pattern to be matched
     */
    void matchStringBidirectionally(Substring& pattern,
                                    std::vector<uint32_t>& nodePathLeft,
                                    std::vector<uint32_t>& nodePathRight);

    /**
     * @brief Static function. Sets the search direction. This way, the
     * correct functions are called.
     *
     * @param d the direction to search in, either FORWARD or BACKWARD
     */
    static void setDirection(Direction d) {
        FMPosSFR::ExtendFMPosAboveK =
            (d == FORWARD) ? &FMPosSFR::extendFMPosAboveKForward
                           : &FMPosSFR::extendFMPosAboveKBackward;
        FMPosSFR::AddCharAboveK = (d == FORWARD)
                                      ? &FMPosSFR::addCharAboveKForward
                                      : &FMPosSFR::addCharAboveKBackward;
    }

    size_t nodePathSize() const {
        return numberOfNodesLeft + numberOfNodesRight;
    }

    /**
     * @brief Operator overloading. Outputs the FMPosSFR object to the
     * outputstream
     *
     * @param os the output stream
     * @param pos the FMPosSFR to print
     * @return std::ostream& - the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const FMPosSFR& pos);
};