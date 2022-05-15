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

#include "fmindex.h"

#include <algorithm> // sorting, minmax
#include <fstream>   // reading in files

using namespace std;

uint k_DBG;

// ============================================================================
// CLASS FMIndex
// ============================================================================

template <class positionClass>
thread_local std::vector<std::vector<FMPosExt<positionClass>>>
    FMIndex<positionClass>::stacks;
template <class positionClass>
thread_local std::vector<BitParallelED> FMIndex<positionClass>::matrices;
template <class positionClass>
thread_local length_t FMIndex<positionClass>::nodeCounter;
template <class positionClass>
thread_local length_t FMIndex<positionClass>::matrixElementCounter;
template <class positionClass>
thread_local length_t FMIndex<positionClass>::positionsInPostProcessingCounter;
template <class positionClass>
thread_local length_t FMIndex<positionClass>::redundantNodePathsCounter;
template <class positionClass>
thread_local Direction FMIndex<positionClass>::dir = BACKWARD;
template <class positionClass>
thread_local ExtraCharPtr<positionClass> FMIndex<positionClass>::extraChar;

// ----------------------------------------------------------------------------
// ROUTINES FOR INITIALIZATION
// ----------------------------------------------------------------------------

template <class positionClass>
void FMIndex<positionClass>::fromFiles(const string& baseFile, bool verbose) {
    if (verbose) {
        cout << "Reading in files with baseFile " << baseFile << endl;

        // read the counts table

        cout << "Reading " << baseFile << ".cct"
             << "\n";
    }

    // read the counts table
    vector<length_t> charCounts(256, 0);
    if (!readArray(baseFile + ".cct", charCounts)) {
        throw runtime_error("Cannot open file: " + baseFile + ".cct");
    }

    // TODO why not in construction process?
    length_t cumCount = 0; // cumulative character counts
    for (size_t i = 0; i < charCounts.size(); i++) {
        if (charCounts[i] == 0)
            continue;
        counts.push_back(cumCount);
        cumCount += charCounts[i];
    }
    sigma = Alphabet<ALPHABET>(charCounts);

    if (verbose) {
        // read the BWT
        cout << "Reading " << baseFile << ".bwt" << endl;
    }
    if (!readText(baseFile + ".bwt", bwt)) {
        throw runtime_error("Cannot open file: " + baseFile + ".bwt");
    }

    textLength = (bwt[bwt.size() - 1] == '\n') ? bwt.size() - 1 : bwt.size();

    if (verbose) {
        cout << "Done reading BWT (size: " << bwt.size() << ")" << endl;

        // read the reverse BWT
        cout << "Reading " << baseFile << ".rev.bwt" << endl;
    }
    if (!readText(baseFile + ".rev.bwt", revbwt)) {
        throw runtime_error("Cannot open file: " + baseFile + ".rev.bwt");
    }
    if (verbose) {
        cout << "Done reading reverse BWT (size: " << revbwt.size() << ")"
             << endl;

        // read the baseFile occurrence table
        cout << "Reading " << baseFile << ".brt" << endl;
    }

    if (!fwdRepr.read(baseFile + ".brt"))
        throw runtime_error("Cannot open file: " + baseFile + ".brt");
    if (verbose) {
        cout << "Done reading baseFile occurrence table" << endl;
    }

    // read the reverse baseFile occurrence table
    if (!revRepr.read(baseFile + ".rev.brt"))
        throw std::runtime_error("Cannot open file: " + baseFile + ".rev.brt");
    if (verbose) {
        cout << "Done reading reverse baseFile occurrence table" << std::endl;
    }
}

template <class positionClass>
void FMIndex<positionClass>::populateTable(bool verbose) {
    if (verbose) {
        cout << "Populating FM-range table with " << wordSize << "-mers...";
    }
    cout.flush();

    Kmer::setWordSize(wordSize);

    table.resize(1 << (2 * wordSize)); // 2 << wordSize is 4^wordSize
    setDirection(FORWARD);

    string word;
    vector<FMPosExt<positionClass>> stack;
    extendFMPos(getCompleteRange(), stack);
    while (!stack.empty()) {
        auto curr = stack.back();
        stack.pop_back();

        word.resize(curr.getRow());
        word[curr.getRow() - 1] = curr.getCharacter();

        if ((length_t)curr.getRow() == wordSize) { // max depth reached
            Kmer k(word);
            table.insert(make_pair(k, curr.getRanges()));

        } else // add extra characters
            extendFMPos(curr.getRanges(), stack, curr.getRow());
    }
    if (verbose) {
        cout << "done." << endl;
    }
}

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING DATA STRUCTURE
// ----------------------------------------------------------------------------

template <class positionClass> void FMIndex<positionClass>::resetCounters() {
    nodeCounter = 0;
    matrixElementCounter = 0;
    positionsInPostProcessingCounter = 0;
    redundantNodePathsCounter = 0;
}

template <class positionClass>
length_t FMIndex<positionClass>::getNodes() const {
    return nodeCounter;
}

template <class positionClass>
length_t FMIndex<positionClass>::getMatrixElements() const {
    return matrixElementCounter;
}

template <class positionClass>
length_t FMIndex<positionClass>::getTotalReported() const {
    return positionsInPostProcessingCounter;
}

template <class positionClass>
length_t FMIndex<positionClass>::getTotalReportedNodePaths() const {
    return redundantNodePathsCounter;
}

template <class positionClass>
length_t FMIndex<positionClass>::findLF(length_t k, bool reversed) const {
    if (reversed) {
        length_t posInAlphabet = sigma.c2i((unsigned char)revbwt[k]);

        return counts[posInAlphabet] + getNumberOfOccRev(posInAlphabet, k);
    }

    length_t posInAlphabet = sigma.c2i((unsigned char)bwt[k]);

    return counts[posInAlphabet] + getNumberOfOcc(posInAlphabet, k);
}

template <class positionClass>
length_t FMIndex<positionClass>::findSA(length_t index) const {
    length_t l = 0;
    while (!sparseSA[index]) {
        index = findLF(index, false);
        l++;
    }
    return sparseSA.get(index) + l;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR EXACT PATTERN MATCHING
// ----------------------------------------------------------------------------
template <class positionClass>
Range FMIndex<positionClass>::matchString(const string& s) {
    // start at the end
    auto it = s.crbegin();

    // find the range for this initial character in the BWT string
    length_t positionInAlphabet = sigma.c2i((unsigned char)*it);

    length_t start = counts[positionInAlphabet];
    length_t end;
    if (positionInAlphabet != sigma.size() - 1) {
        end = counts[positionInAlphabet + 1];
    } else {
        end = bwt.size();
    }
    nodeCounter++;

    // iterate starting from the second character over the string
    for (++it; it != s.crend(); it++) {
        // find number of occurrences of this char before and after and so the
        // new range is found
        positionInAlphabet = sigma.c2i((unsigned char)*it);
        length_t startOfChar = counts[positionInAlphabet];
        start = getNumberOfOcc(positionInAlphabet, start) + startOfChar;
        end = getNumberOfOcc(positionInAlphabet, end) + startOfChar;
        if (start == end) {
            // no matches found
            return Range();
        }
        nodeCounter++;
    }

    // return this range as a pair
    return Range(start, end);
}

template <class positionClass>
vector<length_t> FMIndex<positionClass>::exactMatches(const string& s) {
    // find the range in the suffix array that matches the string
    Range range = matchString(s);

    // declare the return vector
    vector<length_t> positions;
    positions.reserve(range.width());

    // fill in the vector with all values in this range in the suffix array
    for (length_t i = range.getBegin(); i < range.getEnd(); i++) {
        positions.emplace_back(findSA(i));
    }

    // sort the vector and return
    sort(positions.begin(), positions.end());
    return positions;
}

template <class positionClass>
SARangePair
FMIndex<positionClass>::matchStringBidirectionally(const Substring& pattern,
                                                   SARangePair rangesOfPrev) {

    for (length_t i = 0; i < pattern.size(); i++) {

        char c = pattern[i];
        if (!addChar(c, rangesOfPrev)) {
            // rangesOfPrev was made empty
            break;
        }
    }

    return rangesOfPrev;
}

template <class positionClass>
bool FMIndex<positionClass>::addChar(const char& c,
                                     SARangePair& startRange) const {

    int posInAlphabet = sigma.c2i((unsigned char)c);
    if (posInAlphabet > -1) {

        if ((this->*extraChar)(posInAlphabet, startRange, startRange)) {
            // each character that we look at is a new node that is visited
            nodeCounter++;
            return true;
        }
    }
    // the range is now empty

    return false;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

template <class positionClass>
std::vector<TextOccurrence>
FMIndex<positionClass>::approxMatchesNaive(const std::string& pattern,
                                           length_t maxED) {

    std::vector<FMOcc<positionClass>> occurrences =
        approxMatchesNaiveIntermediate(pattern, maxED);

    return mapOccurrencesInSAToOccurrencesInText(occurrences, maxED);
}

template <class positionClass>
std::vector<FMOcc<positionClass>>
FMIndex<positionClass>::approxMatchesNaiveIntermediate(
    const std::string& pattern, length_t maxED) {

    resetCounters();
    vector<FMOcc<positionClass>> occurrences;
    std::vector<uint32_t> leftNodes, rightNodes;

    BandMatrix matrix(pattern.size() + maxED + 1, maxED);

    setDirection(FORWARD);

    std::vector<FMPosExt<positionClass>> stack;
    stack.reserve((pattern.size() + maxED + 1) * (sigma.size() - 1));

    extendFMPos(getCompleteRange(), stack, 0, 0);

    while (!stack.empty()) {
        FMPosExt<positionClass> currentNode = stack.back();
        stack.pop_back();

        updateNodeStackWithNodePath(currentNode, leftNodes, rightNodes);

        int row = currentNode.getDepth();
        length_t first = matrix.getFirstColumn(row);
        length_t last = matrix.getLastColumn(row);
        length_t minimalED = maxED + 1;
        for (length_t j = first; j <= last; j++) {
            matrix.updateMatrix(currentNode.getCharacter() != pattern[j - 1],
                                row, j);
            minimalED = min(minimalED, matrix(row, j));
        }

        if (minimalED > maxED) {
            // backtrack
            continue;
        }

        if (last == (length_t)matrix.getLastColumn()) {
            // full pattern was matched
            if (matrix(row, last) <= maxED) {
                reportMatchEditNaive(currentNode, matrix(row, last),
                                     occurrences, leftNodes, rightNodes);
            }
        }

        extendFMPos(currentNode, stack);
    }

    return occurrences;
}

template <class positionClass>
void FMIndex<positionClass>::setDirection(Direction d) {
    dir = d;
    extraChar = (d == FORWARD) ? &FMIndex::findRangesWithExtraCharForward
                               : &FMIndex::findRangesWithExtraCharBackward;
}

template <class positionClass>
bool FMIndex<positionClass>::findRangesWithExtraCharBackward(
    length_t positionInAlphabet, const SARangePair& rangesOfP,
    SARangePair& rangesOfChild) const {

    // first make the trivial range by  searching cP using B
    Range trivialRange = rangesOfP.getRangeSA();

    // find the new range by using the LF property
    length_t occBefore =
        getNumberOfOcc(positionInAlphabet, trivialRange.getBegin());
    length_t occAfter =
        getNumberOfOcc(positionInAlphabet, trivialRange.getEnd());

    length_t startInAlphabet = counts[positionInAlphabet];
    Range range1 =
        Range(occBefore + startInAlphabet, occAfter + startInAlphabet);

    // then make the less trivial range by counting the sizes of the ranges
    // of (dP) using B

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSARev().getBegin();

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range

    length_t x = getNumberOfCumOcc(positionInAlphabet, trivialRange.getEnd()) -
                 getNumberOfCumOcc(positionInAlphabet, trivialRange.getBegin());

    // make the new range with width equal to that of the trivial range
    Range range2 = Range(s + x, s + x + range1.width());

    rangesOfChild = SARangePair(range1, range2);
    return !rangesOfChild.empty();
}

template <class positionClass>
bool FMIndex<positionClass>::findRangesWithExtraCharForward(
    length_t positionInAlphabet, const SARangePair& rangesOfP,
    SARangePair& childRanges) const {

    // first make the trivial range by searching (Pc)' using B'  if
    // searching forward we need to use B' so we need the reverse range
    Range rangeForTrivial = rangesOfP.getRangeSARev();

    // find the new range by using the LF property
    length_t occBefore =
        getNumberOfOccRev(positionInAlphabet, rangeForTrivial.getBegin());
    length_t occAfter =
        getNumberOfOccRev(positionInAlphabet, rangeForTrivial.getEnd());

    length_t startInAlphabet = counts[positionInAlphabet];
    Range range1 =
        Range(occBefore + startInAlphabet, occAfter + startInAlphabet);

    // then make the less trivial range by counting the size of the range of
    // (Pd)' using B' (forward)

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSA().getBegin();

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range
    Range prevRange = rangesOfP.getRangeSARev();

    length_t x = getNumberOfCumOccRev(positionInAlphabet, prevRange.getEnd()) -
                 getNumberOfCumOccRev(positionInAlphabet, prevRange.getBegin());

    // make the new range
    Range range2 = Range(s + x, s + x + range1.width());

    childRanges = SARangePair(range2, range1);
    return !childRanges.empty();
}

template <class positionClass>
void FMIndex<positionClass>::recApproxMatchEditNaive(
    const Search& s, const FMOcc<positionClass>& startMatch,
    vector<FMOcc<positionClass>>& occ, const vector<Substring>& parts,
    std::vector<uint32_t>& leftNodes, std::vector<uint32_t>& rightNodes,
    const int& idx) {
    const Substring& p = parts[s.getPart(idx)];           // this part
    const length_t& maxED = s.getUpperBound(idx);         // maxED for this part
    const length_t& minED = s.getLowerBound(idx);         // minED for this part
    const length_t& W = maxED - startMatch.getDistance(); // Width of matrix
    const length_t& pSize = p.size();

    BandMatrix matrix = BandMatrix(pSize, W, startMatch.getDistance());

    if (matrix.getLastColumn(0) == (int)pSize && matrix(0, pSize) <= maxED) {
        // an occurrence found by gapping entire part
        if (s.isEnd(idx)) {

            reportMatchEditNaive(startMatch.getPosition(), matrix(0, pSize),
                                 occ, leftNodes, rightNodes);
        } else {
            // go to the next index
            recApproxMatchEditNaive(
                s,
                FMOcc<positionClass>(startMatch.getPosition(),
                                     matrix(0, pSize)),
                occ, parts, leftNodes, rightNodes, idx + 1);
            // set direction correct again
            setDirection(dir);
        }
    }

    auto& stack = stacks[idx];                  // stack for this partition
    const Direction& dir = s.getDirection(idx); // direction
    setDirection(dir);

    positionClass copy = startMatch.getPosition();
    copy.setDepth(0);
    extendFMPos(copy, stack);

    while (!stack.empty()) {
        auto currentNode = stack.back();
        stack.pop_back();

        updateNodeStackWithNodePath(currentNode, leftNodes, rightNodes);

        const int& row = currentNode.getRow();
        const length_t firstCol = matrix.getFirstColumn(row);
        const length_t lastCol = matrix.getLastColumn(row);

        length_t minimalEDOfRow = matrix(row, firstCol - 1);
        for (length_t col = firstCol; col <= lastCol; col++) {
            length_t filledIn = matrix.updateMatrix(
                currentNode.getCharacter() != p[col - 1], row, col);
            matrixElementCounter++;
            minimalEDOfRow = std::min(minimalEDOfRow, filledIn);
        }

        if (minimalEDOfRow > maxED) {
            // backtracking
            continue;
        }

        if (lastCol == pSize && matrix(row, pSize) <= maxED &&
            matrix(row, pSize) >= minED) {
            length_t originalDepth = currentNode.getDepth();
            currentNode.setDepth(startMatch.getDepth() +
                                 currentNode.getDepth());
            if (s.isEnd(idx)) {
                reportMatchEditNaive(currentNode, matrix(row, pSize), occ,
                                     leftNodes, rightNodes);
            } else {
                // go deeper in search
                recApproxMatchEditNaive(
                    s, FMOcc<positionClass>(currentNode, matrix(row, pSize)),
                    occ, parts, leftNodes, rightNodes, idx + 1);
                // set direction correct again
                setDirection(dir);
            }
            currentNode.setDepth(originalDepth);
        }
        if (firstCol == pSize) {
            // final cell of matrix filled in => backtracking
            continue;
        }

        // extend to the children
        extendFMPos(currentNode, stack);
    }
}

template <class positionClass>
void FMIndex<positionClass>::recApproxMatchEditOptimized(
    const Search& s, const FMOcc<positionClass>& startMatch,
    vector<FMOcc<positionClass>>& occ, const vector<Substring>& parts,
    vector<uint32_t>& leftNodes, vector<uint32_t>& rightNodes, const int& idx,
    const vector<FMPosExt<positionClass>>& descPrevDir,
    const vector<uint>& initPrevDir,
    const vector<FMPosExt<positionClass>>& descNotPrevDir,
    const vector<uint>& initNotPrevDir) {

    // shortcut Variables
    const Substring& p = parts[s.getPart(idx)];      // this part
    const length_t& maxED = s.getUpperBound(idx);    // maxED for this part
    const Direction& dir = s.getDirection(idx);      // direction
    const bool& dSwitch = s.getDirectionSwitch(idx); // has direction switched?
    auto& stack = stacks[idx];                       // stack for this partition
    size_t matrixIdx = s.getPart(idx) + (dir == BACKWARD) * s.getNumParts();
    BitParallelED& bpED = matrices[matrixIdx]; // matrix for this partition

    // get the correct initED and descendants based on switch
    const vector<uint>& initEds = dSwitch ? initNotPrevDir : initPrevDir;
    const vector<FMPosExt<positionClass>>& descendants =
        dSwitch ? descNotPrevDir : descPrevDir;
    const vector<uint>& initOther = dSwitch ? initPrevDir : initNotPrevDir;
    const vector<FMPosExt<positionClass>>& descOther =
        dSwitch ? descPrevDir : descNotPrevDir;

    setDirection(dir);

    // calculate necessary increase for first column of bandmatrix
    vector<uint> initED;
    if (initEds.empty()) {
        initED = vector<uint>(1, startMatch.getDistance());
    } else {
        uint prevED = (dSwitch ? *min_element(initEds.begin(), initEds.end())
                               : initEds[0]);
        uint increase = startMatch.getDistance() - prevED;
        initED = vector<uint>(initEds.size());
        for (size_t i = 0; i < initED.size(); i++) {
            initED[i] = initEds[i] + increase;
        }
    }

    stack.reserve((p.size() + maxED + 1) * ALPHABET);

    // encode the sequence of this partition in the matrix if this has not been
    // done before
    if (!bpED.sequenceSet())
        bpED.setSequence(p);

    // initialize bit-parallel matrix
    bpED.initializeMatrix(maxED, initED);

    // initialize cluster
    Cluster<positionClass> clus(bpED.getSizeOfFinalColumn(), maxED,
                                startMatch.getDepth(), startMatch.getShift());

    if (bpED.inFinalColumn(0)) {
        // the first row is part of the final column of the banded matrix
        // Update the first cell of the cluster with the startmatch and the
        // value found at the psizeth column of the initialization row the
        // character does not matter as the first cell of a cluster is never a
        // descendant of any of the other cells in the cluster, so this cell
        // will not be reused for the next part of the pattern
        positionClass copy = startMatch.getPosition();
        copy.setDepth(0);
        clus.setValue(0, FMPosExt<positionClass>((char)0, copy),
                      bpED(0, p.size()));
    }

    if (!descendants.empty()) {
        // fill in the matrix for the descendants

        length_t maxRow = bpED.getNumberOfRows() - 1;

        for (length_t i = 0;
             i < descendants.size() && descendants[i].getDepth() <= maxRow;
             i++) {

            if (branchAndBound(
                    clus, descendants[i], s, idx, parts, occ, leftNodes,
                    rightNodes, initOther, descOther,
                    {descendants.begin() + i + 1, descendants.end()})) {
                return;
            }
        }
        if (descendants.back().getDepth() == maxRow) {
            // minimal ed exceeded no more options to get lower ed from
            // first column, or no more rows to possibly check
            return;
        }

        if (dSwitch) {
            // after a switch the ranges of final descendant should be
            // updated
            positionClass copy = startMatch.getPosition();
            copy.setDepth(descendants.back().getDepth());
            // push children of final descendant
            extendFMPos(copy, stack);
        } else {
            // push children of final descendant
            extendFMPos(descendants.back(), stack);
        }
    } else { // get the initial nodes to check
        positionClass startMatchCopy = startMatch.getPosition();
        startMatchCopy.setDepth(0);
        extendFMPos(startMatchCopy, stack);
    }

    while (!stack.empty()) {

        FMPosExt<positionClass> currentNode = stack.back();
        stack.pop_back();

        updateNodeStackWithNodePath(currentNode, leftNodes, rightNodes);

        if (branchAndBound(clus, currentNode, s, idx, parts, occ, leftNodes,
                           rightNodes, initOther, descOther)) {

            continue;
        }

        // continue the search for children of this node
        extendFMPos(currentNode, stack);
    }
}

template <class positionClass>
const bool FMIndex<positionClass>::separationIsNext(positionClass pos) const {
    for (int i = 0; i < numberOfSeparationCharacters; i++) {
        SARangePair pairForNewChar;
        if ((this->*extraChar)(i, pos.getRanges(), pairForNewChar)) {
            return true;
        }
    }
    return false;
}

template <class positionClass>
bool FMIndex<positionClass>::branchAndBound(
    Cluster<positionClass>& clus, const FMPosExt<positionClass>& currentNode,
    const Search& s, const int& idx, const vector<Substring>& parts,
    vector<FMOcc<positionClass>>& occ, vector<uint32_t>& leftNodes,
    vector<uint32_t>& rightNodes, const vector<uint>& initOther,
    const vector<FMPosExt<positionClass>>& descOther,
    const vector<FMPosExt<positionClass>> remainingDesc) {

    // get the appropriate matrix
    size_t matrixIdx = s.getPart(idx) + (dir == BACKWARD) * s.getNumParts();
    BitParallelED& bpED = matrices[matrixIdx];

    const Substring& p = parts[s.getPart(idx)]; // The current partition
    const length_t& pSize = p.size(); // The size of the current partition

    // compute, in a bit-parallel manner, a single row of the ED matrix
    const length_t row = currentNode.getDepth();
    bool validED = bpED.computeRow(row, currentNode.getCharacter());

    // check if we have reached the final column of the matrix
    const length_t lastCol = bpED.getLastColumn(row);

    if (lastCol == pSize) {
        // update the cluster
        size_t clusIdx = clus.size() + row - bpED.getNumberOfRows();
        clus.setValue(clusIdx, currentNode, bpED(row, lastCol));

        // During strain-free matching, we need to go deeper into the matrix,
        // even if this part could only be reached with vertical gaps. This is
        // necessary because all possible node paths need to be reported for a
        // complete filtering process later on.
        if (!validED || (bpED.onlyVerticalGapsLeft(row) && !strainFree) ||
            (clusIdx + 1 == clus.size() && strainFree)) {
            // no need to further explore this branch for this part -> go to
            // next part
            goDeeper(clus, idx + 1, s, parts, occ, leftNodes, rightNodes,
                     s.getLowerBound(idx), descOther, initOther, remainingDesc);
            return true;
        }

        // TODO reduce redundancy here!
        if (s.isEdge(idx) && separationIsNext(currentNode)) {
            // no need to further explore this branch for this part -> go to
            // next part
            goDeeper(clus, idx + 1, s, parts, occ, leftNodes, rightNodes,
                     s.getLowerBound(idx), descOther, initOther, remainingDesc);
            // See if it is also followed by other characters, so no return yet
        }
    }

    return !validED;
}

template <class positionClass>
void FMIndex<positionClass>::goDeeper(
    Cluster<positionClass>& cluster, const length_t& nextIdx, const Search& s,
    const vector<Substring>& parts, vector<FMOcc<positionClass>>& occ,
    vector<uint32_t>& leftNodes, vector<uint32_t>& rightNodes,
    const length_t& lowerBound,
    const vector<FMPosExt<positionClass>>& descOtherD,
    const vector<uint>& initOtherD,
    const vector<FMPosExt<positionClass>>& remainingDesc) {

    bool isEdge = s.isEdge(nextIdx - 1);

    if (isEdge) {
        // if this is final piece report highest minimum (to get shortest
        // match)
        if (nextIdx == parts.size()) {
            vector<FMOcc<positionClass>> matches;
            cluster.reportCentersAtEnd(matches, leftNodes, rightNodes);
            for (const auto& match : matches) {
                if (match.isValid() && match.getDistance() >= (int)lowerBound) {
                    occ.emplace_back(match);
                }
            }
        } else {
            std::vector<FMOcc<positionClass>> matches;
            cluster.reportNonFinalEdgePiece(this->dir, matches);
            Direction originalDir = this->dir;
            // go deeper in search for every reported match
            for (auto match : matches) {
                if (match.isValid() && match.getDistance() >= (int)lowerBound) {
                    // make a copy the left and right nodes as we might need
                    // these originals for backtracking in the CURRENT part
                    auto left = leftNodes, right = rightNodes;
                    recApproxMatchEditOptimized(s, match, occ, parts, left,
                                                right, nextIdx, {}, {},
                                                descOtherD, initOtherD);
                }
            }
            // set direction back again
            setDirection(originalDir);
        }

        return;
    }

    // one of the later stages will return to this point, so keep track of
    // the descendants and eds at this branch
    vector<FMPosExt<positionClass>> descendants;
    vector<uint> initEds;

    FMOcc<positionClass> newMatch =
        cluster.getClusterCentra(lowerBound, descendants, initEds);
    if (!newMatch.isValid()) {
        // no centre above lowerbound found
        return;
    }

    // add the remaining descendants
    descendants.insert(descendants.end(), remainingDesc.begin(),
                       remainingDesc.end());

    // reset the depth of all the descendants
    for (size_t i = 0; i < descendants.size(); i++) {
        descendants[i].setDepth(i + 1);
    }
    uint maxEDNext = s.getUpperBound(nextIdx);

    // remove trailing initEds that are higher than maxEDNext
    while (initEds.back() > maxEDNext) {
        initEds.pop_back();
    }

    // is the next direction equal to this direction?
    bool switchAfter = s.getDirectionSwitch(nextIdx);

    if (switchAfter) {
        // switching direction as this is not the end of a search direction,
        // this means we'll get back here, thus range of newmatch should be
        // deepest point in branch
        if (!descendants.empty()) {
            length_t originalDepth = newMatch.getPosition().getDepth();
            newMatch.setPosition(descendants.back());
            newMatch.setDepth(originalDepth);

            //  edit distance for search in other direction should be lowest
            //  value possible
            newMatch.setDistance(*min_element(initEds.begin(), initEds.end()));
        }

        Direction originalDir = this->dir;

        recApproxMatchEditOptimized(s, newMatch, occ, parts, leftNodes,
                                    rightNodes, nextIdx, descendants, initEds,
                                    descOtherD, initOtherD);

        // set direction back again
        setDirection(originalDir);
    } else {
        // go deeper on next piece
        recApproxMatchEditOptimized(s, newMatch, occ, parts, leftNodes,
                                    rightNodes, nextIdx, descendants, initEds,
                                    descOtherD, initOtherD);
    }
}

// TODO: make this function compatible with strain-free matching
template <class positionClass>
void FMIndex<positionClass>::recApproxMatchHamming(
    const Search& s, const FMOcc<positionClass>& startMatch,
    std::vector<FMOcc<positionClass>>& occ, const std::vector<Substring>& parts,
    const int& idx) {

    // shortcut variables
    const Substring& p = parts[s.getPart(idx)]; // the current part
    const int& pSize = p.size();                // the size of the current part
    const Direction& d = s.getDirection(idx);   // direction of current part
    const int& maxED = s.getUpperBound(idx);    // upperbound of current part
    const int& minED = s.getLowerBound(idx);    // lowerbound of currentpart
    setDirection(d);

    // create vector
    std::vector<int> vec(p.size() + 1, 0);
    // set root element of vector to the distance of the startmatch
    vec[0] = startMatch.getDistance();
    // get stack for current part
    auto& stack = stacks[idx];

    extendFMPos(startMatch.getRanges(), stack);

    while (!stack.empty()) {
        FMPosExt<positionClass> node = stack.back();
        stack.pop_back();

        // update the vector
        int row = node.getRow();
        vec[row] = vec[row - 1] + (node.getCharacter() != p[row - 1]);
        if (vec[row] > maxED) {
            // backtrack
            continue;
        }

        if (row == pSize) {
            // end of part
            if (vec[row] >= minED) {
                // valid occurrence
                FMOcc<positionClass> match = FMOcc<positionClass>(
                    node.getRanges(), vec[row], startMatch.getDepth() + pSize);
                if (s.isEnd(idx)) {
                    // end of search
                    occ.emplace_back(match);
                } else {
                    // continue search
                    recApproxMatchHamming(s, match, occ, parts, idx + 1);
                }
            }
            continue;
        }
        extendFMPos(node, stack);
    }
}

template <class positionClass>
bool FMIndex<positionClass>::extendFMPosIntermediary(
    const SARangePair& parentRanges, vector<FMPosExt<positionClass>>& stack,
    int row, length_t i, int trueDepth) {

    SARangePair pairForNewChar;

    // check if this character occurs in the specified range
    if ((this->*extraChar)(i, parentRanges, pairForNewChar)) {
        // push this range and character for the next iteration
        stack.emplace_back(sigma.i2c(i), pairForNewChar, row + 1,
                           trueDepth + 1);

        nodeCounter++;
        return true;
    }
    return false;
}

template <class positionClass>
void FMIndex<positionClass>::extendFMPos(const SARangePair& parentRanges,
                                         vector<FMPosExt<positionClass>>& stack,
                                         int row, int trueDepth) {

    // iterate over the entire alphabet
    for (length_t i = 2; i < sigma.size(); i++) {

        extendFMPosIntermediary(parentRanges, stack, row, i);
    }
}

template <class positionClass>
void FMIndex<positionClass>::extendFMPos(
    const positionClass& pos, vector<FMPosExt<positionClass>>& stack) {
    extendFMPos(pos.getRanges(), stack, pos.getDepth());
}

// ----------------------------------------------------------------------------
// POST-PROCESSING ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

template <class positionClass>
vector<TextOccurrence> FMIndex<positionClass>::convertToMatchesInText(
    const FMOcc<positionClass>& saMatch) {

    vector<TextOccurrence> textMatches;
    textMatches.reserve(saMatch.getRanges().width());

    for (length_t i = saMatch.getRanges().getRangeSA().getBegin();
         i < saMatch.getRanges().getRangeSA().getEnd(); i++) {
        // find the startPosition in the text by looking at the SA
        length_t startPos = findSA(i) + saMatch.getShift();

        // cap startPos at textLength
        startPos = startPos % textLength;

        length_t endPos = startPos + saMatch.getDepth();

        textMatches.emplace_back(Range(startPos, endPos),
                                 saMatch.getDistance());
    }
    return textMatches;
}

template <class positionClass>
vector<TextOccurrence>
FMIndex<positionClass>::mapOccurrencesInSAToOccurrencesInText(
    vector<FMOcc<positionClass>>& occ, const int& maxED) {

    sort(occ.begin(), occ.end());
    occ.erase(unique(occ.begin(), occ.end()), occ.end());
    vector<TextOccurrence> occurrencesInText;
    occurrencesInText.reserve(1000 * maxED);

    if (occ.size() == 0) {
        return {};
    }
    if (occ.size() == 1) {
        // all occ are distinct
        positionsInPostProcessingCounter = occ[0].getRanges().width();
        auto m = convertToMatchesInText(occ[0]);
        sort(m.begin(), m.end());
        for (auto& occ : m) {
            occ.generateOutput();
        }
        return m;
    }

    // more than 1 reported occurrence, could be redundant
    for (const auto& it : occ) {

        const Range& range = it.getRanges().getRangeSA();
        positionsInPostProcessingCounter += range.width();

        auto matchesInTextToCheck = convertToMatchesInText(it);

        occurrencesInText.insert(occurrencesInText.end(),
                                 matchesInTextToCheck.begin(),
                                 matchesInTextToCheck.end());
    }

    sort(occurrencesInText.begin(), occurrencesInText.end());

    std::vector<TextOccurrence> nonRedundantOcc;
    nonRedundantOcc.reserve(occurrencesInText.size());

    length_t maxDiff = 2 * maxED;
    length_t prevBegin = numeric_limits<length_t>::max();
    int prevED = maxED + 1;
    length_t prevDepth = numeric_limits<length_t>::max();

    for (const auto& o : occurrencesInText) {
        auto diff = abs_diff<length_t>(o.getRange().getBegin(), prevBegin);
        if (diff == 0) {
            continue;
        }
        if (diff <= maxDiff) {
            // check if this later occurrence is better than the previous
            // one
            if (o.getDistance() > prevED) {
                continue;
            }
            if (o.getDistance() == prevED &&
                o.getRange().width() >= prevDepth) {
                continue;
            }

            // prev was worse so pop_back
            nonRedundantOcc.pop_back();
        }

        prevBegin = o.getRange().getBegin();
        prevED = o.getDistance();
        prevDepth = o.getRange().width();

        nonRedundantOcc.emplace_back(o);
    }
    for (TextOccurrence& occ : nonRedundantOcc) {
        occ.generateOutput();
    }

    return nonRedundantOcc;
}

template class FMIndex<FMPos>;
template class FMIndex<FMPosSFR>;