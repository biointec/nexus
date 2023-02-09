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

#ifndef FMINDEX_H
#define FMINDEX_H

#include "alphabet.h"
#include "bandmatrix.h"
#include "buildIndexAuxiliary.h"
#include "bwtrepr.h"
#include "cluster.h"
#include "encodedtext.h"
#include "search.h"
#include "suffixarray.h"
#include "textoccurrence.h"
#include "tkmer.h"

#include <google/sparse_hash_map>
#include <string>
#include <vector>

// ============================================================================
// (TYPE) DEFINITIONS AND PROTOTYPES
// ============================================================================

typedef uint64_t length_t;

// ============================================================================
// CLASS FMIndex
// ============================================================================

template <class positionClass>
using ExtraCharPtr = bool (FMIndex<positionClass>::*)(length_t,
                                                      const SARangePair&,
                                                      SARangePair&) const;

template <class positionClass> class FMIndex {
    friend class FMPos;

  protected:
    // info about the text

    //   basefile of the reference text
    const std::string baseFile;
    // The length of the reference text
    length_t textLength;
    // The reference text
    EncodedText<ALPHABET> text;

    // The number of separation characters
    int numberOfSeparationCharacters = 1;

    // the alphabet
    Alphabet<ALPHABET> sigma; // the alphabet

    // bidirectional FM-index data structures

    // the bwt string of the reference genome
    EncodedText<ALPHABET> bwt;
    // the bwt string of the reverse reference genome
    EncodedText<ALPHABET> revbwt;
    // the counts array of the reference genome
    std::vector<length_t> counts;
    // the (sparse) suffix array of the reference genome
    SparseSuffixArray sparseSA;
    // the baseFile occurrences table
    BWTRepr<ALPHABET> fwdRepr;
    // the baseFile occurrences of the rev BWT
    BWTRepr<ALPHABET> revRepr;

    // performance counters

    // Count the number of nodes visited
    thread_local static length_t nodeCounter;
    // Count the number of matrix elements visited
    thread_local static length_t matrixElementCounter;
    // Count the number of positions that are post processed
    thread_local static length_t positionsInPostProcessingCounter;
    // Count the number of reported node paths
    thread_local static length_t redundantNodePathsCounter;

    // direction variables

    // the direction of the index
    thread_local static Direction dir;
    // pointer to extra char method (for direction)
    thread_local static ExtraCharPtr<positionClass> extraChar;

    // stacks for search schemes: stacks of nodes for the different partitions
    thread_local static std::vector<std::vector<FMPosExt<positionClass>>>
        stacks;
    thread_local static std::vector<BitParallelED>
        matrices; // alignment matrices for the different partitions

    // sparse hash info

    // the size the mers to be stored in a table
    static const size_t wordSize = 4;
    // hashtable that contains all wordSize-mers
    google::sparse_hash_map<Kmer, SARangePair, KmerHash> table;

    // Indicates whether the matching is strain-free
    bool strainFree = false;

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * @brief Construct a new FMIndex object. This private constructor is only
     * used in the building process. It uses a stub parameter to distinguish
     * itself from other constructors. Except for storing the baseFile string,
     * nothing happens in this constructor.
     *
     * @param stub Stub parameter to create a different constructor for the
     * building process. Set this value to 0.
     * @param baseFile the baseFile of the files that will be read in
     * @param sa_sparse sparseness factor of suffix array. It is assumed this is
     * a power of two
     * @param verbose if true the steps will be written to cout
     */
    FMIndex(int stub, const std::string& baseFile, int sa_sparse = 1,
            bool verbose = true)
        : baseFile(baseFile) {
    }

    /**
     * @brief Private helper function that reads in all the necessary files to
     * construct an FM-index that has already been built.
     *
     * @param baseFile the baseFile of the files that will be read in
     * @param verbose if true the steps will be written to cout
     * @param strainFree bool indicating whether strainfree matching is required
     */
    void fromFiles(const std::string& baseFile, bool verbose,
                   bool strainFree = false);

    /**
     * @brief Read a binary file and stores content in array (e.g. suffix array)
     *
     * @param filename File name
     * @param array Suffix array (contents will be overwritten)
     * @return true if successful
     * @return false otherwise
     */
    static bool readArray(const std::string& filename,
                          std::vector<length_t>& array) {
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs)
            return false;

        ifs.seekg(0, std::ios::end);
        array.resize(ifs.tellg() / sizeof(length_t));
        ifs.seekg(0, std::ios::beg);
        ifs.read((char*)&array[0], array.size() * sizeof(length_t));

        return true;
    }

    /**
     * @brief Read a text file (e.g. input text, BWT, ...)
     *
     * @param filename File name
     * @param buf Buffer (contents will be overwritten)
     * @return true if successful
     * @return false otherwise
     */
    static bool readText(const std::string& filename, std::string& buf) {
        std::ifstream ifs(filename);
        if (!ifs)
            return false;

        ifs.seekg(0, std::ios::end);
        buf.resize(ifs.tellg());
        ifs.seekg(0, std::ios::beg);
        ifs.read((char*)&buf[0], buf.size());

        return true;
    }

    /**
     * @brief Populate the hash table
     *
     * @param verbose if true, the steps are written to cout
     */
    void populateTable(bool verbose);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * @brief Finds the LF mapping of the character at index k in the (reverse)
     * bwt string
     *
     * @param k the index to find the LF mapping of
     * @param reversed If true, the LF property is applied on the reverse BWT.
     * If false, it is applied on the regular BWT.
     * @return length_t - the row that is the LF mapping of k over the (reverse)
     * BWT. It is so that the entry in the (reverse) suffix array of this return
     * value is one less than the entry in the (reverse) suffix array at index k
     */
    length_t findLF(length_t k, bool reversed) const;

    /**
     * @brief Find the entry in the suffix array of this index. This is computed
     * from the sparse suffix array and the bwt
     *
     * @param index the index to find the entry in the SA of
     * @return length_t - the entry in the SA of the index
     */
    length_t findSA(length_t index) const;

    /**
     * @brief Function that returns the number of occurrences in the bwt before
     * an index of the symbol at symbolindex in the alphabet
     *
     * @param symbolIndex the index in the alphabet to get the number of
     * occurrences of
     * @param index the index whose entry for symbol in the occurrences table is
     * asked
     * @return length_t - the number of occurrences of the symbol before index
     * in the bwt
     */
    length_t getNumberOfOcc(length_t symbolIndex, length_t index) const {
        return fwdRepr.occ(symbolIndex, index);
    }

    /**
     * @brief Function that returns the number of occurrences in the reverse bwt
     * before an index of the symbol at symbolindex in the alphabet
     *
     * @param symbolIndex the index in the alphabet to get the number of
     * occurrences of
     * @param index the index whose entry for symbol in the reverse occurrences
     * table is asked
     * @return length_t - the number of occurrences of the symbol before index
     * in the reverse bwt
     */
    length_t getNumberOfOccRev(length_t symbolIndex, length_t index) const {
        return revRepr.occ(symbolIndex, index);
    }

    /**
     * @brief Function that returns the number of occurrences in the bwt before
     * an index of all symbols smaller than the symbol at symbolindex in the
     * alphabet
     *
     * @param symbolIndex the index in the alphabet to get the cumulative number
     * of occurrences of
     * @param index the index whose entry for symbol in the prefixoccurrences
     * table is asked
     * @return length_t - the number of occurrences of symbols smaller than
     * symbol at symbolindex before index index in the bwt
     */
    length_t getNumberOfCumOcc(length_t symbolIndex, length_t index) const {
        return fwdRepr.cumOcc(symbolIndex, index);
    }

    /**
     * @brief Function that returns the number of occurrences in the reverse bwt
     * before an index of all symbols smaller than the symbol at symbolindex in
     * the alphabet
     *
     * @param symbolIndex the index in the alphabet to get the cumulative number
     * of occurrences of
     * @param index the index whose entry for symbol in the reverse
     * prefixoccurrences table is asked
     * @return length_t - length_t - the number of occurrences of symbols
     * smaller than symbol at symbolindex before index index in the reverse bwt
     */
    length_t getNumberOfCumOccRev(length_t symbolIndex, length_t index) const {
        return revRepr.cumOcc(symbolIndex, index);
    }

    // ----------------------------------------------------------------------------
    // HELP ROUTINES FOR EXACT MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Private helper function for exactMatches. This finds the range in
     * the suffix array that corresponds to matches in the reference genome
     *
     * @param s the string to be matched in the reference genome
     * @return Range - a Range containing the start and end values of the range
     * ([start, end[)
     */
    Range matchString(const std::string& s);

    // ----------------------------------------------------------------------------
    // HELPER ROUTINES FOR APPROXIMATE MATCHING (ITERATIVELY)
    // ----------------------------------------------------------------------------

    /**
     * @brief Goes deeper in a search if a valid approximate match is found in
     * the cluster
     *
     * @param cluster the cluster to search for a valid approximate
     * @param nextp the  idx of next part to research
     * @param s the search
     * @param parts the parts of the pattern
     * @param occ the vector with all occurrences of the entire pattern, if the
     * current partition is the final partition of the search then the match (if
     * one found) will be pushed onto this vector
     * @param leftNodes the left nodepath (from origin to leftmost point)
     * @param rightNodes the right node path (from origin to rightmost point)
     * @param lowerBound the lowerbound for this partition
     * @param descendantsOtherD the descendants of the other direction, defaults
     * to empty vector
     * @param initEdsOtherD the initialization eds of the other direction,
     * defaults to empty vector
     * @param remainingDesc the remaining descendants on the current branch,
     * that are already created but aren't checked yet and need to be checked
     * for the next part, defaults to an empty vector
     */
    void
    goDeeper(Cluster<positionClass>& cluster, const length_t& nextp,
             const Search& s, const std::vector<Substring>& parts,
             std::vector<FMOcc<positionClass>>& occ,
             std::vector<uint32_t>& leftNodes,
             std::vector<uint32_t>& rightNodes, const length_t& lowerBound,
             const std::vector<FMPosExt<positionClass>>& descendantsOtherD = {},
             const std::vector<uint>& initEdsOtherD = {},
             const std::vector<FMPosExt<positionClass>>& remainingDesc = {});

    /**
     * @brief Check if one of the children of a position is a separation
     * character.
     *
     * @param pos The position of which the children need to be checked
     * @return true if one of the children of the position is a separation
     * character
     * @return false otherwise
     */
    virtual const bool separationIsNext(positionClass pos) const;

    /**
     * @brief Helper function for the approximate matching. This function fills
     * in the matrix for the current node at the current row and goes deeper for
     * the next part is that is necessary
     *
     * @param matrix the matrix to fill in
     * @param clus the cluster corresponding to the final column of the matrix
     * @param currentNode the node for which the matrix is filled in
     * @param s the current search
     * @param idx the idx of the current part
     * @param parts the parts of the pattern
     * @param occ vector to push occurrences to
     * @param leftNodes the left nodepath (from origin to leftmost point)
     * @param rightNodes the right node path (from origin to rightmost point)
     * @param initOther eds of descendants in other direction
     * @param descOther descendants in other direction
     * @param remainingDesc the remaining descendants on the current branch,
     * that are already created but aren't checked yet and might need to be
     * checked for the next part, defaults to an empty vector
     * @return true if the search can backtrack
     * @return false false if the search can continue along this branch for the
     * current part
     */
    bool branchAndBound(
        Cluster<positionClass>& clus,
        const FMPosExt<positionClass>& currentNode, const Search& s,
        const int& idx, const std::vector<Substring>& parts,
        std::vector<FMOcc<positionClass>>& occ,
        std::vector<uint32_t>& leftNodes, std::vector<uint32_t>& rightNodes,
        const std::vector<uint>& initOther,
        const std::vector<FMPosExt<positionClass>>& descOther,
        const std::vector<FMPosExt<positionClass>> remainingDesc = {});

    /**
     * Dummy function for template argument
     * @param pos the final position
     * @param leftNodes the left nodepath (from origin to leftmost point)
     * @param rightNodes the right node path (from origin to rightmost point)
     */
    void updateNodeStackWithNodePath(const FMPosExt<FMPos>& pos,
                                     std::vector<uint32_t>& leftNodes,
                                     std::vector<uint32_t>& rightNodes) {
        // this method does nothing as the path does not need to be updated as
        // there is no path
    }

    /**
     * Updates the node path stacks according to the current position. The new
     * path is always either a subpath of the current path or the current path
     * with one extra node. The entire path is reversed(leftNodes) + rightNodes
     * @param pos the current position
     * @param leftNodes the left nodepath (from origin to leftmost point)
     * @param rightNodes the right node path (from origin to rightmost point)
     */
    void updateNodeStackWithNodePath(const FMPosExt<FMPosSFR>& pos,
                                     std::vector<uint32_t>& leftNodes,
                                     std::vector<uint32_t>& rightNodes) {
        pos.updateNodeStackWithNodePath(leftNodes, rightNodes);
    }

    /**
     * @brief Pushes the child that was found by extending node with ranges
     * equal to ranges with the ith character onto the stack
     *
     * @param parentRanges the ranges to get the child of
     * @param stack the stack to push the child on
     * @param row the row of the parentNode
     * @param i the rank of the character to add
     * @param trueDepth the true depth of the total current match (defaults
     * to -1)
     * @return true if a new extension was pushed onto the stack
     * @return false otherwise
     */
    bool extendFMPosIntermediary(const SARangePair& parentRanges,
                                 std::vector<FMPosExt<positionClass>>& stack,
                                 int row, length_t i, int trueDepth = -1) const;

    /**
     * @brief Pushes all the children corresponding to the node with ranges
     * equal to ranges onto the stack
     *
     * @param ranges the ranges to get the children of
     * @param stack the stack to push the children on
     * @param row the row of the parentNode (defaults to 0)
     * @param trueDepth the true depth of the total current match (defaults to
     * -1)
     */
    virtual void extendFMPos(const SARangePair& ranges,
                             std::vector<FMPosExt<positionClass>>& stack,
                             int row = 0, int trueDepth = -1) const;

    /**
     * @brief Pushes all the children corresponding to the this position onto
     * the stack
     *
     * @param pos the position to get the children of
     * @param stack the stack to push the children on
     */
    virtual void extendFMPos(const positionClass& pos,
                             std::vector<FMPosExt<positionClass>>& stack) const;

    /**
     * @brief Converts a match in the suffix array to matches in the text
     *
     * @param matchInSA the match that will be converted
     * @return std::vector<TextOccurrence> - a vector with the corresponding
     * text occurrences
     */
    std::vector<TextOccurrence>
    convertToMatchesInText(const FMOcc<positionClass>& matchInSA);

    /**
     * @brief Add a new occurrence to the list of matches
     *
     * @param pos the position object of the new occurrence
     * @param distance the edit distance of the new occurrence
     * @param occ the list of matches to add to
     * @param leftNodes empty list of nodes
     * @param rightNodes empty list of nodes
     */
    void reportMatchEditNaive(FMPos pos, int distance,
                              std::vector<FMOcc<FMPos>>& occ,
                              std::vector<uint32_t>& leftNodes,
                              std::vector<uint32_t>& rightNodes) {
        occ.emplace_back(pos, distance);
    }

    /**
     * @brief Set the correct node path pointer for the position object and add
     * the new occurrence to the list of matches
     *
     * @param pos the position object of the new occurrence. Node path pointer
     * will be added
     * @param distance the edit distance of the new occurrence
     * @param occ the list of matches to add to
     * @param leftNodes left side of the node path, reverse
     * @param rightNodes right side of the node path
     */
    void reportMatchEditNaive(FMPosSFR pos, int distance,
                              std::vector<FMOcc<FMPosSFR>>& occ,
                              std::vector<uint32_t>& leftNodes,
                              std::vector<uint32_t>& rightNodes) {

        // find nodePath at this particular node
        // start by making a copy of left
        auto path = leftNodes;

        // constrict left to  appropriate number of nodes
        path.resize(pos.getNumberOfNodesLeft());

        // reverse left
        std::reverse(path.begin(), path.end());
        // add right to left
        path.insert(path.end(), rightNodes.begin(),
                    rightNodes.begin() + pos.getNumberOfNodesRight());

        pos.setNodePath(path);

        occ.emplace_back(pos, distance);
    }

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * @brief Construct a new FMIndex object that has already been built
     *
     * @param baseFile baseFile of the files that contain the info
     * @param sa_sparse sparseness factor of suffix array. It is assumed this is
     * a power of two
     * @param verbose will write to cout
     * @param strainFree bool indicating whether strain free matching is
     * required
     */
    FMIndex(const std::string& baseFile, int sa_sparse = 1, bool verbose = true,
            bool strainFree = false)
        : baseFile(baseFile), sparseSA(baseFile, sa_sparse) {
        // read in files
        fromFiles(baseFile, verbose, strainFree);

        // populate table
        populateTable(verbose);
    }

    /**
     * @brief Get the complete range of this index
     *
     * @return SARangePair - an SARangePair with both ranges the complete range
     * of the index
     */
    SARangePair getCompleteRange() const {
        return SARangePair(Range(0, bwt.size()), Range(0, bwt.size()));
    }

    /**
     * @brief Get an empty SARangePair object. Specifically, all bounds are set
     * to 0.
     *
     * @return SARangePair - an empty SARangePair object
     */
    SARangePair getEmptyRange() const {
        return SARangePair();
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING THE DATASTRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * @brief Get a reference to the original text
     *
     * @return const EncodedText<ALPHABET>& - a reference to the original text
     */
    const EncodedText<ALPHABET>& getText() {
        if (text.empty()) {
            if (!text.read(baseFile + ".compressed.txt")) {
                throw std::runtime_error("Problem reading: " + baseFile +
                                         ".txt");
            }
        }
        return text;
    }

    /**
     * @brief Reset the counters
     *
     */
    void resetCounters();

    /**
     * @brief Get the node counter
     *
     * @return length_t - the node counter
     */
    length_t getNodes() const;

    /**
     * @brief Get the matrix element counter
     *
     * @return length_t - the matrix element counter
     */
    length_t getMatrixElements() const;

    /**
     * @brief Get the counter for the number of positions in post processing
     *
     * @return length_t - the counter for the number of positions in post
     * processing
     */
    length_t getTotalReported() const;

    /**
     * @brief Get the counter for the number of reported node paths
     *
     * @return length_t - the number of reported node paths
     */
    length_t getTotalReportedNodePaths() const;

    /**
     * @brief Get the wordsize of the mers stored in the table
     *
     * @return int - the wordsize of the mers stored in the table
     */
    int getWordSize() const {
        return wordSize;
    }

    /**
     * @brief Get the ranges of a single character in this index
     *
     * @param c the character we want to find the ranges for
     * @return SARangePair - the ranges of character c in this index
     */
    SARangePair getRangeOfSingleChar(char c) const {
        unsigned int i = sigma.c2i(c);
        if (i < sigma.size() - 1) {
            return SARangePair(Range(counts[i], counts[i + 1]),
                               Range(counts[i], counts[i + 1]));
        }
        return SARangePair(Range(counts[i], bwt.size()),
                           Range(counts[i], bwt.size()));
    }

    /**
     * @brief Looks up the SARangePair corresponding to k-mer p in the
     * hashtable. Assumes p is of size wordSize.
     *
     * @param p the substring to find the ranges of
     * @return SARangePair - the ranges corresponding to substring p, if no pair
     * can be found returns empty ranges
     */
    SARangePair lookUpInKmerTable(const Substring& p) const {
        Kmer k(p.tostring());

        auto it = table.find(k);
        if (it != table.end()) {
            return it->second;
        }

        return SARangePair();
    }

    /**
     * @brief Get the text length
     *
     * @return length_t - text length
     */
    const length_t& getTextLength() const {
        return textLength;
    }

    /**
     * @brief Get the counts
     *
     * @return std::vector<length_t>& - counts
     */
    const std::vector<length_t>& getCounts() const {
        return counts;
    }

    /**
     * @brief Get the alphabet
     *
     * @return const Alphabet<ALPHABET> - the alphabet
     */
    const Alphabet<ALPHABET> getSigma() const {
        return sigma;
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR EXACT MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Calculate the positions in the reference genome where exact
     * matches to the argument string start.
     *
     * @param s the string to match in the reference genome
     * @return std::vector<length_t> - a sorted vector containing the start
     * positions of all exact substring matches of s in the reference sequence
     */
    std::vector<length_t> exactMatches(const std::string& s);

    /**
     * @brief Matches a string exactly and return the ranges in the SA and
     * reverse SA
     *
     * @param pattern the string to match
     * @return SARangePair - the pair of ranges of this pattern
     */
    SARangePair matchStringBidirectionally(const Substring& pattern) {
        return matchStringBidirectionally(pattern, getCompleteRange());
    }

    /**
     * @brief Matches a string exactly starting form startRange
     *
     * @param pattern the string to match
     * @param startRange the range to search in
     * @return SARangePair - the resulting pair of ranges
     */
    SARangePair matchStringBidirectionally(const Substring& pattern,
                                           SARangePair startRange);

    /**
     * @brief Add one character and update the range. If the character cannot be
     * added the range will be set to an empty range.
     *
     * @param c the character to be added (in the current direction of the
     * index)
     * @param range the range to extend, will be overwritten
     * @return true if the resulting range is not empty
     * @return false otherwise
     */
    bool addChar(const char& c, SARangePair& range) const;

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Resizes to the required number of stacks and reserves space on
     * each stack, so that each stack can match the entire pattern
     *
     * @param number the number of stacks required
     * @param size the size of the pattern
     */
    void reserveStacks(const length_t number, const int size) {
        stacks.resize(number);
        length_t stackSize = size * sigma.size();
        for (auto& stack : stacks) {
            stack.reserve(stackSize);
        }
    }

    /**
     * Reset the in-text matrices to be empty matrices
     * @param number the number of partitions
     */
    void resetMatrices(const length_t number) {
        matrices.resize(2 * number);

        for (auto& matrix : matrices) {
            matrix.reset();
        }
    }

    // ----------------------------------------------------------------------------
    // HELP ROUTINES FOR APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Find the ranges of cP using the principle explained in the paper
     * of Lam
     *
     * @param positionInAlphabet the position in the alphabet of the character
     * that is added in the front
     * @param rangesOfP the ranges of pattern P
     * @param childRanges the ranges cP, this will be overwritten
     * @return true if the new ranges are not empty
     * @return false otherwise
     */
    bool findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                         const SARangePair& rangesOfP,
                                         SARangePair& childRanges) const;

    /**
     * @brief Find the ranges of Pc using the principle explained in the paper
     * of Lam
     *
     * @param positionInAlphabet the position in the alphabet of the character
     * that is added in the back
     * @param rangesOfP the ranges of pattern P
     * @param childRanges the ranges Pc, this will be overwritten
     * @return true if the new ranges are not empty
     * @return false otherwise
     */
    bool findRangesWithExtraCharForward(length_t positionInAlphabet,
                                        const SARangePair& rangesOfP,
                                        SARangePair& childRanges) const;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Match the pattern approximately in a naive way. All matches are at
     * most a certain edit distance away from the pattern. This function also
     * maps the occurrences in the SA to occurrences in the original text.
     *
     * @param pattern the pattern to match
     * @param maxED the maximum edit distance
     * @return std::vector<TextOccurrence> - a vector with matches which contain
     * a range (the range of the text that matched) and the edit distance this
     * substring is away from the pattern
     */
    std::vector<TextOccurrence> approxMatchesNaive(const std::string& pattern,
                                                   length_t maxED);

    /**
     * @brief Private helper function for the naive approximate pattern matching
     * method
     *
     * @param pattern the pattern to match
     * @param maxED the maximum edit distance
     * @return std::vector<FMOcc<positionClass>> - a vector with occurrences in
     * the FM-index or the SA
     */
    std::vector<FMOcc<positionClass>>
    approxMatchesNaiveIntermediate(const std::string& pattern, length_t maxED);

    /**
     * @brief Sets the search direction of the fm-index
     *
     * @param d the direction to search in, either FORWARD or BACKWARD
     */
    virtual void setDirection(Direction d);

    /**
     * @brief Matches a search recursively with a depth first approach (each
     * branch of the tree is fully examined until the backtracking condition is
     * met) using hamming distance metric
     *
     * @param s the search to follow
     * @param startMatch the approximate match found for all previous parts of
     * the search
     * @param occ a vector with matches of the complete search, if a such a
     * match is found it is pushed upon this vector
     * @param parts the parts of the pattern that needs to be matched
     * @param idx the index of the partition to match, defaults to 1 as an exact
     * search for the zeroth part is assumed
     */
    void recApproxMatchHamming(const Search& s,
                               const FMOcc<positionClass>& startMatch,
                               std::vector<FMOcc<positionClass>>& occ,
                               const std::vector<Substring>& parts,
                               const int& idx = 1);

    /**
     * @brief Matches a search recursively with a depth first approach (each
     * branch of the tree is fully examined until the backtracking condition is
     * met) using the edit distance metric. This function uses all optimizations
     * for eliminating redundancy in the edit distance metric
     *
     * @param search the search to follow
     * @param startMatch the approximate match found for all previous partitions
     * of the search
     * @param occ a vector with matches of the complete search, if such a match
     * is found is a pushed upon this vector
     * @param parts the parts of the pattern that needs to be matched
     * @param leftNodes the left nodepath (from origin to leftmost point)
     * @param rightNodes the right node path (from origin to rightmost point)
     * @param idx the index of the partition to match, defaults to 1 as an exact
     * search for the zeroth partition is assumed
     * @param descPrevDir the descendants of the previous direction, defaults to
     * empty vector
     * @param initPrevDir the initialization eds of the previous direction,
     * defaults to empty vector
     * @param descNotPrevDir the descendants of the other direction, defaults to
     * empty vector
     * @param initNotPrevDir the initialization eds of the other direction,
     * defaults to empty vector
     */
    void recApproxMatchEditOptimized(
        const Search& search, const FMOcc<positionClass>& startMatch,
        std::vector<FMOcc<positionClass>>& occ,
        const std::vector<Substring>& parts, std::vector<uint32_t>& leftNodes,
        std::vector<uint32_t>& rightNodes, const int& idx = 1,
        const std::vector<FMPosExt<positionClass>>& descPrevDir =
            std::vector<FMPosExt<positionClass>>(),
        const std::vector<uint>& initPrevDir = std::vector<uint>(),
        const std::vector<FMPosExt<positionClass>>& descNotPrevDir =
            std::vector<FMPosExt<positionClass>>(),
        const std::vector<uint>& initNotPrevDir = std::vector<uint>());

    /**
     * @brief Matches a search recursively with a depth first approach (each
     * branch of the tree is fully examined until the backtracking condition is
     * met) using the edit distance metric. This function does not use any
     * optimizations for eliminating redundancy in the edit distance metric. It
     * simply matches the current part starting from startrange and each node
     * found that has an edit distance between the lower and upperbound is used
     * to start a search for the next part
     *
     * @param s the search to follow
     * @param startMatch the approximate match found for all previous partitions
     * of the search
     * @param occ a vector with matches of the complete search, if such a match
     * is found is a pushed upon this vector
     * @param parts the parts of the pattern that needs to be matched
     * @param idx the index of the partition to match, defaults to 1 as an exact
     * @param leftNodes the left nodepath (from origin to leftmost point)
     * @param rightNodes the right node path (from origin to rightmost point)
     * search for the zeroth partition is assumed
     */
    void recApproxMatchEditNaive(const Search& s,
                                 const FMOcc<positionClass>& startMatch,
                                 std::vector<FMOcc<positionClass>>& occ,
                                 const std::vector<Substring>& parts,
                                 std::vector<uint32_t>& leftNodes,
                                 std::vector<uint32_t>& rightNodes,
                                 const int& idx);

    // ----------------------------------------------------------------------------
    // POST-PROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief This function maps matches in the SA to matches in the text. It
     * takes the ranges of the matches and together with the depth this is
     * matched to a range in the text (this new range has a width of depth). The
     * edit distance is also mapped to this new range in the text. This function
     * also filters out the redundant matches using the maximal allowed edit
     * distance
     *
     * @param occurrences the vector with all approximate matches and their
     * ranges in the SA and revSA
     * @param maxED the maximal allowed edit distance
     * @return std::vector<TextOccurrence> - a vector of matches in the text
     * containing a range and the edit distance
     */
    std::vector<TextOccurrence> mapOccurrencesInSAToOccurrencesInText(
        std::vector<FMOcc<positionClass>>& occurrences, const int& maxED);
};

#endif
