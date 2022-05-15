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

#include "fmpos.h"

#include <vector>

// compute |a-b| in a safe manner
template <typename T> T abs_diff(T a, T b) {
    return a > b ? a - b : b - a;
}

// ============================================================================
// CLASS FMOcc
// ============================================================================

/**
 * An occurrence in the bidirectional FM-index
 */
template <class positionClass> class FMOcc {

  protected:
    // the position in the FMIndex of this approximate match
    positionClass pos;
    // the distance (hamming or edit)
    int distance;
    // A right-sift to the corresponding positions in the text
    length_t shift;

  public:
    /**
     * @brief Construct a new FMOcc object. Default constructor: distance = 0
     * and shift = 0.
     *
     */
    FMOcc() : pos(), distance(0), shift(0) {
    }

    /**
     * @brief Construct a new FMOcc object or a bidirectional approximate match
     * in the suffix array.
     *
     * @param ranges the ranges of this approximate match (range in SA and in
     * reverse SA)
     * @param distance the (edit or hamming) distance of this approximate match
     * @param depth the depth (=length) of this approximate match
     * @param shift the right shift to the corresponding positions in the text,
     * defaults to zero
     */
    FMOcc(SARangePair ranges, int distance, int depth, int shift = 0)
        : pos(ranges, depth), distance(distance),
          shift(shift) { // TODO should true depth be set here as well?
    }

    /**
     * @brief Construct a new FMOcc object or a bidirectional approximate match
     * in the suffix array
     *
     * @param pos the position in the FMIndex of this approximate match
     * @param distance the (edit or hamming) distance of this approximate match
     * @param shift The right shift to the corresponding positions in the text,
     * defaults to zero
     */
    FMOcc(positionClass pos, int distance, int shift = 0)
        : pos(pos), distance(distance), shift(shift) {
    }

    /**
     * @brief Get the ranges of this approximate match (range in SA and in
     * reverse SA)
     *
     * @return const SARangePair& - the ranges of this approximate match (range
     * in SA and in reverse SA)
     */
    const SARangePair& getRanges() const {
        return pos.getRanges();
    }

    /**
     * @brief Get the (edit or hamming) distance of this approximate match
     *
     * @return const int& - the (edit or hamming) distance of this approximate
     * match
     */
    const int& getDistance() const {
        return distance;
    }

    /**
     * @brief Get the depth (=length) of this approximate match
     *
     * @return const length_t& - the depth (=length) of this approximate match
     */
    const length_t& getDepth() const {
        return pos.getDepth();
    }

    /**
     * @brief Get the right shift to the corresponding positions in the text
     *
     * @return const length_t& - the right shift to the corresponding positions
     * in the text
     */
    const length_t& getShift() const {
        return shift;
    }

    /**
     * @brief Get the position in the FMIndex of this approximate match
     *
     * @return const positionClass& - the position in the FMIndex of this
     * approximate match
     */
    const positionClass& getPosition() const {
        return pos;
    }

    /**
     * @brief Set the ranges of this approximate match (range in SA and in
     * reverse SA)
     *
     * @param ranges the ranges of this approximate match (range in SA and in
     * reverse SA)
     */
    void setRanges(SARangePair ranges) {
        pos.setRanges(ranges);
    }

    /**
     * @brief Set the (edit or hamming) distance of this approximate match
     *
     * @param distance the (edit or hamming) distance of this approximate match
     */
    void setDistance(int distance) {
        this->distance = distance;
    }

    /**
     * @brief Set the depth (=length) of this approximate match
     *
     * @param depth the depth (=length) of this approximate match
     */
    void setDepth(length_t depth) {
        pos.setDepth(depth);
    }

    /**
     * @brief Set the true depth of this approximate match
     *
     * @param trueDepth the true depth of this approximate match
     */
    void setTrueDepth(length_t trueDepth) {
        pos.setTrueDepth(trueDepth);
    }

    /**
     * @brief Set the position in the FMIndex of this approximate match
     *
     * @param pos the position in the FMIndex of this approximate match
     */
    void setPosition(positionClass pos) {
        this->pos = pos;
    }

    /**
     * @brief Check if the position in the FMIndex of this approximate match is
     * valid
     *
     * @return true if the position is valid
     * @return false otherwise
     */
    bool isValid() const {
        return pos.isValid();
    }

    /**
     * @brief Operator overloading to sort FMOcc. The specific sorting technique
     * depends on the class type of the position object: FMPos or FMPosSFR.
     *
     * @param rhs the FMOcc to compare to this
     * @return true if this is smaller than rhs
     * @return false otherwise
     */
    bool operator<(const FMOcc<positionClass>& rhs) const {

        return pos.compare(rhs, distance, shift);
    }

    /**
     * @brief Operator overloading. Two FMocc are equal if their position,
     * distance and shift are all equal.
     *
     * @param rhs
     * @return true if this is equal to rhs
     * @return false otherwise
     */
    bool operator==(const FMOcc<positionClass>& rhs) {
        return getPosition() == rhs.getPosition() &&
               distance == rhs.getDistance() && getShift() == rhs.getShift();
    }

    /**
     * @brief Operator overloading. Outputs the FMOcc object to the outputstream
     *
     * @param os the output stream
     * @param biocc the FMOcc to print
     * @return std::ostream& - the output stream
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const FMOcc<FMPos>& biocc);

    /**
     * @brief Operator overloading. Outputs the FMOcc object to the outputstream
     *
     * @param os the output stream
     * @param biocc the FMOcc to print
     * @return std::ostream& - the output stream
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const FMOcc<FMPosSFR>& biocc);
};

// ============================================================================
// CLASS FMOccSFI
// ============================================================================

template <class positionClass> class FMOccSFI : public FMOcc<positionClass> {

    // This is class is used to contain occurrences in the FM-index along
    // with their location in the graph. This class is only used during
    // strain-fixed matching.

  private:
    // The node path in the graph that corresponds to this FMOcc
    std::vector<uint32_t> nodepath;

  public:
    /**
     * @brief Construct a new FMOccSFI object
     *
     * @param occ The FMOcc object on which this object will be built
     * @param nodepath The node path in the graph to be added
     */
    FMOccSFI(FMOcc<positionClass> occ, const std::vector<uint32_t>& nodepath)
        : FMOcc<positionClass>(occ), nodepath(nodepath) {
    }

    /**
     * @brief Get the node path in the graph
     *
     * @return const std::vector<uint32_t>& - the node path in the graph
     */
    const std::vector<uint32_t>& getNodePath() const {
        return nodepath;
    }
};

// ============================================================================
// CLASS FMOccSFR
// ============================================================================

// Declaration necessary for typedefs below
class FMOccSFR;

// Create function pointer types
typedef void (*replacesPtr)(FMOccSFR* thisOcc, FMOccSFR* replacedOcc);
typedef bool (*setReplacementPtr)(FMOccSFR* thisOcc, FMOccSFR* replacementOcc);
typedef void (*reverseNodePathPtr)(FMOccSFR& thisOcc);
typedef void (*reportPtr)(FMOccSFR* thisOcc);

class FMOccSFR : public FMOcc<FMPosSFR> {
    // These occurrences are only used during strain-free matching
  private:
    // ----------------------------------------------------------------------------
    // ATTRIBUTES FOR COMPLETE FILTERING
    // ----------------------------------------------------------------------------

    // A bool indicating whether the occurrence must not be removed during
    // filtering, because it serves as a replacement for an occurrence that was
    // removed earlier.
    bool toBeKept = false;
    // A bool indicating that the occurrence might have to be removed during
    // filtering, because it there is another occurrence better than this one.
    bool toBeRemoved = false;

    // ----------------------------------------------------------------------------
    // ATTRIBUTES FOR LINEAR FILTERING
    // ----------------------------------------------------------------------------

    // Bool indicating whether the occurrence is reported when the regular node
    // path is considered.
    bool reportedRegular = false;
    // Bool indicating whether the occurrence is reported when the reverse node
    // path is considered.
    bool reportedReverse = false;

    // Stores the depth of the shortest occurrence that is replaced by this
    // occurrence when the regular node path is considered.
    length_t depthOfShortestReplacedRegular;
    // Stores the depth of the longest occurrence that is replaced by this
    // occurrence when the regular node path is considered.
    length_t depthOfLongestReplacedRegular;
    // Stores the depth of the shortest occurrence that is replaced by this
    // occurrence when the reverse node path is considered.
    length_t depthOfShortestReplacedReverse;
    // Stores the depth of the longest occurrence that is replaced by this
    // occurrence when the reverse node path is considered.
    length_t depthOfLongestReplacedReverse;

    // Stores the edit distance of the best replacement for this occurrence when
    // the regular node path is considered.
    int EDofReplacementRegular;
    // Stores the edit distance of the best replacement for this occurrence when
    // the reverse node path is considered.
    int EDofReplacementReverse;

    // Stores the depth of the best replacement for this occurrence when the
    // regular node path is considered.
    length_t depthOfReplacementRegular;
    // Stores the depth of the best replacement for this occurrence when the
    // reverse node path is considered.
    length_t depthOfReplacementReverse;

  public:
    /**
     * @brief Empty constructor. Used to make dummy objects.
     *
     */
    FMOccSFR() : FMOcc() {
    }

    /**
     * @brief Construct a new FMOccSFR object or a strain-free approximate
     * matche in the pan-genome compressed de Bruijn graph
     *
     * @param pos the position in the pan-genome compressed de Bruijn graph of
     * this approximate match
     * @param distance the (edit or hamming) distance of this approximate match
     * @param shift The right shift to the corresponding positions in the graph,
     * defaults to zero
     */
    FMOccSFR(FMPosSFR pos, int distance, int shift = 0)
        : FMOcc(pos, distance, shift) {
        depthOfShortestReplacedRegular = getDepth();
        depthOfLongestReplacedRegular = depthOfShortestReplacedRegular;
        depthOfShortestReplacedReverse = depthOfShortestReplacedRegular;
        depthOfLongestReplacedReverse = depthOfShortestReplacedRegular;
        EDofReplacementRegular = distance;
        EDofReplacementReverse = distance;
        depthOfReplacementRegular = depthOfShortestReplacedRegular;
        depthOfReplacementReverse = depthOfShortestReplacedRegular;
    }

    /**
     * @brief Construct a new FMOccSFR object based on a FMOcc<FMPosSFR> object.
     *
     * @param oldPos The FMOcc<FMPosSFR> object from which the new FMOccSFR
     * object must be built.
     */
    FMOccSFR(const FMOcc<FMPosSFR>& oldPos) : FMOcc(oldPos) {
        depthOfShortestReplacedRegular = getDepth();
        depthOfLongestReplacedRegular = depthOfShortestReplacedRegular;
        depthOfShortestReplacedReverse = depthOfShortestReplacedRegular;
        depthOfLongestReplacedReverse = depthOfShortestReplacedRegular;
        EDofReplacementRegular = distance;
        EDofReplacementReverse = distance;
        depthOfReplacementRegular = depthOfShortestReplacedRegular;
        depthOfReplacementReverse = depthOfShortestReplacedRegular;
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR COMPLETE FILTERING
    // ----------------------------------------------------------------------------

    /**
     * @brief Check whether this occurrence must not be removed during
     * filtering, because it serves as a replacement for an occurrence that was
     * removed earlier.
     *
     * @return true if this occurrence must be kept
     * @return false otherwise
     */
    const bool isToBeKept() const {
        return toBeKept;
    }

    /**
     * @brief Check whether this occurrence might have to be removed during
     * filtering, because it there is another occurrence better than this one.
     *
     * @return true if this occurrence must be kept
     * @return false otherwise
     */
    const bool isToBeRemoved() const {
        return toBeRemoved;
    }

    /**
     * @brief Mark that this occurrence must not be removed during filtering,
     * because it serves as a replacement for an occurrence that was removed
     * earlier.
     *
     */
    void setToBeKept() {
        toBeKept = true;
    }

    /**
     * @brief Mark that this occurrence might have to be removed during
     * filtering, because it there is another occurrence better than this one.
     *
     */
    void setToBeRemoved() {
        toBeRemoved = true;
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR LINEAR FILTERING
    // ----------------------------------------------------------------------------

    /**
     * @brief Static function that sets the direction in which the node path is
     * observed for filtering: regular or reverse.
     *
     * @param reverse True for reverse, false for regular
     */
    void static setDirection(bool reverse) {
        // Set the function pointers to the correct values
        replaces =
            reverse ? &FMOccSFR::replacesReverse : &FMOccSFR::replacesRegular;
        setReplacement = reverse ? &FMOccSFR::setReplacementReverse
                                 : &FMOccSFR::setReplacementRegular;
        reverseNodePath = reverse ? &FMOccSFR::reverseNodePathReverse
                                  : &FMOccSFR::reverseNodePathRegular;
        report = reverse ? &FMOccSFR::reportReverse : &FMOccSFR::reportRegular;
    }

    /**
     * @brief Static function. Object relpacement replaces thisOcc. Therefore,
     * the replacements attributes of thisOcc must be updated based on object
     * replacement. This is only done if object replacement is effectively
     * better than any other replacement that was found up until now (in both
     * directions).
     *
     * @param thisOcc the occurrence to be updated
     * @param replacement the replacing occurrence
     * @return true
     * @return false
     */
    static setReplacementPtr setReplacement;

    /**
     * @brief Static function. Object relpacement replaces thisOcc. Therefore,
     * the replacements attributes of thisOcc must be updated based on object
     * replacement. This is only done if object replacement is effectively
     * better than any other replacement that was found up until now (in both
     * directions). This function is for the reverse direction.
     *
     * @param thisOcc the occurrence to be updated
     * @param replacement the replacing occurrence
     * @return true
     * @return false
     */
    static bool setReplacementReverse(FMOccSFR* thisOcc,
                                      FMOccSFR* replacement) {
        // Store some variables so the data of the pointers does not need to be
        // reached multiple times
        length_t replacementDepth = replacement->getDepth();
        int replacementDistance = replacement->getDistance();
        length_t depthOfReplacementReverse = thisOcc->depthOfReplacementReverse;
        int EDofReplacementReverse = thisOcc->EDofReplacementReverse;
        length_t thisDepth = thisOcc->getDepth();
        length_t depthOfReplacementRegular = thisOcc->depthOfReplacementRegular;
        int EDofReplacementRegular = thisOcc->EDofReplacementRegular;

        // First check that the depth of the replacement is different of that of
        // this object. If this is not the case, replacement and thisOcc are the
        // same which is not of interest. Next, check that the replacement is
        // effectively better than any replacement that was found earlier for
        // this direction. This is done based on the distance and depth values.
        if (replacementDepth != thisDepth &&
            ((replacementDepth >= depthOfReplacementReverse &&
              replacementDistance < EDofReplacementReverse) ||
             (replacementDepth < depthOfReplacementReverse &&
              replacementDistance <= EDofReplacementReverse))) {
            // We have effectively found a better replacement in this direction

            // Set the new depth and edit distance values for the replacement of
            // thisOcc
            thisOcc->depthOfReplacementReverse = replacementDepth;
            thisOcc->EDofReplacementReverse = replacementDistance;

            if (replacementDepth > thisDepth) {
                // If the replacement is longer than thisOcc, it also takes over
                // the shorter occurrences that were previously replaced by
                // thisOcc. Therefore, the attribute that stores the depth of
                // the shortest occurrence that was replaced by thisOcc is
                // reset.
                thisOcc->depthOfShortestReplacedReverse = thisDepth;
            }

            // Now check if the replacement that was found is also better in the
            // other direction. Only if this is the case, the attributes of the
            // replacement occurrence should be updated.
            if ((replacementDepth >= depthOfReplacementRegular &&
                 replacementDistance < EDofReplacementRegular) ||
                (replacementDepth < depthOfReplacementRegular &&
                 replacementDistance <= EDofReplacementRegular)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Static function. Object relpacement replaces thisOcc. Therefore,
     * the replacements attributes of thisOcc must be updated based on object
     * replacement. This is only done if object replacement is effectively
     * better than any other replacement that was found up until now (in both
     * directions). This function is for the regular direction.
     *
     * @param thisOcc the occurrence to be updated
     * @param replacement the replacing occurrence
     * @return true
     * @return false
     */
    static bool setReplacementRegular(FMOccSFR* thisOcc,
                                      FMOccSFR* replacement) {
        // Store some variables so the data of the pointers does not need to be
        // reached multiple times
        length_t replacementDepth = replacement->getDepth();
        int replacementDistance = replacement->getDistance();
        length_t depthOfReplacementReverse = thisOcc->depthOfReplacementReverse;
        int EDofReplacementReverse = thisOcc->EDofReplacementReverse;
        length_t thisDepth = thisOcc->getDepth();
        length_t depthOfReplacementRegular = thisOcc->depthOfReplacementRegular;
        int EDofReplacementRegular = thisOcc->EDofReplacementRegular;

        // First check that the depth of the replacement is different of that of
        // this object. If this is not the case, replacement and thisOcc are the
        // same which is not of interest. Next, check that the replacement is
        // effectively better than any replacement that was found earlier for
        // this direction. This is done based on the distance and depth values.
        if (replacementDepth != thisDepth &&
            ((replacementDepth >= depthOfReplacementRegular &&
              replacementDistance < EDofReplacementRegular) ||
             (replacementDepth < depthOfReplacementRegular &&
              replacementDistance <= EDofReplacementRegular))) {
            // We have effectively found a better replacement in this direction

            // Set the new depth and edit distance values for the replacement of
            // thisOcc
            thisOcc->depthOfReplacementRegular = replacementDepth;
            thisOcc->EDofReplacementRegular = replacementDistance;

            if (replacementDepth > thisDepth) {
                // If the replacement is longer than thisOcc, it also takes over
                // the shorter occurrences that were previously replaced by
                // thisOcc. Therefore, the attribute that stores the depth of
                // the shortest occurrence that was replaced by thisOcc is
                // reset.
                thisOcc->depthOfShortestReplacedRegular = thisDepth;
            }

            // Now check if the replacement that was found is also better in the
            // other direction. Only if this is the case, the attributes of the
            // replacement occurrence should be updated.
            if ((replacementDepth >= depthOfReplacementReverse &&
                 replacementDistance < EDofReplacementReverse) ||
                (replacementDepth < depthOfReplacementReverse &&
                 replacementDistance <= EDofReplacementReverse)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Static function. This function is called if thisOcc replaces
     * replacedOcc. It makes sure that the corresponding attributes in thisOcc
     * are updated.
     *
     * @param thisOcc The replacing occurrence. This occurrence should be
     * updated.
     * @param replacedOcc The occurrence that is being replaced.
     */
    static replacesPtr replaces;

    /**
     * @brief Static function. This function is called if thisOcc replaces
     * replacedOcc. It makes sure that the corresponding attributes in thisOcc
     * are updated. This function is for the reverse direction.
     *
     * @param thisOcc The replacing occurrence. This occurrence should be
     * updated.
     * @param replacedOcc The occurrence that is being replaced.
     */
    static void replacesReverse(FMOccSFR* thisOcc, FMOccSFR* replacedOcc) {
        // Store variables to avoid repeated dereferencing
        length_t replacedDepth = replacedOcc->getDepth();

        if (thisOcc->getDepth() < replacedDepth) {
            // The replaced occurrence is longer than this occurrence
            thisOcc->depthOfLongestReplacedReverse =
                std::max(thisOcc->depthOfLongestReplacedReverse, replacedDepth);
        } else {
            // The replaced occurrence is shorter than this occurrence
            thisOcc->depthOfShortestReplacedReverse = std::min(
                thisOcc->depthOfShortestReplacedReverse, replacedDepth);
        }
    }

    /**
     * @brief Static function. This function is called if thisOcc replaces
     * replacedOcc. It makes sure that the corresponding attributes in thisOcc
     * are updated. This function is for the regular direction.
     *
     * @param thisOcc The replacing occurrence. This occurrence should be
     * updated.
     * @param replacedOcc The occurrence that is being replaced.
     */
    static void replacesRegular(FMOccSFR* thisOcc, FMOccSFR* replacedOcc) {
        // Store variables to avoid repeated dereferencing
        length_t replacedDepth = replacedOcc->getDepth();
        if (thisOcc->getDepth() < replacedDepth) {
            // The replaced occurrence is longer than this occurrence
            thisOcc->depthOfLongestReplacedRegular =
                std::max(thisOcc->depthOfLongestReplacedRegular, replacedDepth);
        } else {
            // The replaced occurrence is shorter than this occurrence
            thisOcc->depthOfShortestReplacedRegular = std::min(
                thisOcc->depthOfShortestReplacedRegular, replacedDepth);
        }
    }

    /**
     * @brief Static function. Reverse the node path and reset the necessary
     * attributes.
     *
     * @param thisOcc
     */
    static reverseNodePathPtr reverseNodePath;

    /**
     * @brief Static function. Reverse the node path and reset the necessary
     * attributes. This function is for the reverse direction.
     *
     * @param thisOcc
     */
    static void reverseNodePathReverse(FMOccSFR& thisOcc) {
        // Reverse the node path
        thisOcc.pos.reverseNodePath();
        // Set distanceFromLeftEnd to distanceFromRightEnd and vice versa
        thisOcc.pos.swapDistances();

        // Store variables to avoid repeated dereferencing
        length_t depth = thisOcc.getDepth();

        // Reset the reverse attributes
        thisOcc.reportedReverse = false;
        thisOcc.depthOfShortestReplacedReverse = depth;
        thisOcc.depthOfLongestReplacedReverse = depth;
        thisOcc.EDofReplacementReverse = thisOcc.getDistance();
        thisOcc.depthOfReplacementReverse = depth;
    }

    /**
     * @brief Static function. Reverse the node path and reset the necessary
     * attributes. This function is for the regular direction.
     *
     * @param thisOcc
     */
    static void reverseNodePathRegular(FMOccSFR& thisOcc) {
        // Reverse the node path
        thisOcc.pos.reverseNodePath();
        // Set distanceFromLeftEnd to distanceFromRightEnd and vice versa
        thisOcc.pos.swapDistances();

        // Store variables to avoid repeated dereferencing
        length_t depth = thisOcc.getDepth();

        // Reset the reverse attributes
        thisOcc.reportedRegular = false;
        thisOcc.depthOfShortestReplacedRegular = depth;
        thisOcc.depthOfLongestReplacedRegular = depth;
        thisOcc.EDofReplacementRegular = thisOcc.getDistance();
        thisOcc.depthOfReplacementRegular = depth;
    }

    /**
     * @brief Static function. This function sets the reported bool of the
     * object parameter to true.
     *
     * @param thisOcc The occurrence for which the reported attribute must be
     * set to true
     */
    static reportPtr report;

    /**
     * @brief Static function. This function sets the reported bool of the
     * object parameter to true. This function is for the reverse direction.
     *
     * @param thisOcc The occurrence for which the reported attribute must be
     * set to true
     */
    static void reportReverse(FMOccSFR* thisOcc) {
        thisOcc->reportedReverse = true;
    }

    /**
     * @brief Static function. This function sets the reported bool of the
     * object parameter to true. This function is for the regular direction.
     *
     * @param thisOcc The occurrence for which the reported attribute must be
     * set to true
     */
    static void reportRegular(FMOccSFR* thisOcc) {
        thisOcc->reportedRegular = true;
    }

    /**
     * @brief Check if two occurrences are located close enough in the
     * pan-genome graph such that they can possibly replace each other. This is
     * based on the maximum allowed edit distance.
     *
     * @param thisOcc Occurrence 1 to check proximity for
     * @param replacementOcc Occurrence 2 to check proximity for
     * @param maxED the maximum edit distance allowed
     * @return true if the occurrences are close enough to replace each other
     * @return false otherwise
     */
    static bool checkProximity(FMOccSFR* thisOcc, FMOccSFR* replacementOcc,
                               length_t maxED) {
        length_t maxDiff = 2 * maxED;
        auto diff = abs_diff<length_t>(
            thisOcc->getPosition().getDistanceFromLeftEnd(),
            replacementOcc->getPosition().getDistanceFromLeftEnd());
        auto temp = diff <= maxDiff;
        return temp;
    }

    /**
     * @brief This function checks if this occurrence is non-redundant according
     * to the linear filtering system.
     *
     * @return true if this occurrence is considered non-redundant
     * @return false otherwise
     */
    const bool isNonRedundant() const {
        // Store the depth of this occurrence
        length_t depth = getDepth();

        if ((reportedRegular && reportedReverse)) {
            // If the occurrence is reported in both directions, it the minimum
            // for a prefix branch in both directions. In this case, the
            // occurrence is non-redundant.
            return true;
        } else if ((reportedRegular &&
                    ((depthOfShortestReplacedRegular >= depth &&
                      depth > depthOfReplacementReverse) ||
                     (depthOfLongestReplacedRegular <= depth &&
                      depth < depthOfReplacementReverse))) ||
                   (reportedReverse &&
                    ((depthOfShortestReplacedReverse >= depth &&
                      depth > depthOfReplacementRegular) ||
                     (depthOfLongestReplacedReverse <= depth &&
                      depth < depthOfReplacementRegular)))) {
            // If the occurrence is only reported in one direction, then it is a
            // minimum in a prefix branch in this direction, but has a better
            // replacement in the other direction. If however, the replacement
            // is shorter than this occurrence and all replaced occurrences in
            // the other direction are longer than this occurrence, than the
            // replacement from the first direction also replaces the longer
            // occurrences in the other direction. For this reason, this
            // occurrence is redundant after all. The same goes for a longer
            // replacement and shorter replaced occurrences.
            return false;
        } else if ((reportedRegular && depthOfShortestReplacedRegular !=
                                           depthOfLongestReplacedRegular) ||
                   (reportedReverse && depthOfShortestReplacedReverse !=
                                           depthOfLongestReplacedReverse)) {
            // Here, the occurrence is only reported in one direction, then it
            // is a minimum in a prefix branch in this direction, but has a
            // better replacement in the other direction. The case that was
            // explained above does not apply however. In this case, there is no
            // efficient way of knowing whether this occurrence is redundant or
            // not, so it is marked as not-redundant.
            return true;
        } else {
            // Leftover cases are marked as redundant.
            return false;
        }
    }
};