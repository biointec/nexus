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

#ifndef SEARCHSTRATEGY_H
#define SEARCHSTRATEGY_H

#include "fmindexDBG.h"

#define Pattern std::vector<int>

// An enum for partion strategy
enum PartitionStrategy { UNIFORM, STATIC, DYNAMIC };
// An enum for which distance metric to use
enum DistanceMetric { HAMMING, EDITNAIVE, EDITOPTIMIZED };

// ============================================================================
// CLASS SEARCHSTRATEGY
// ============================================================================

// This is an abstract class. Every derived class should be able to create
// searches for a given value of k. This abstract base class handles the
// partitioning (either with values provided in the derived class or default
// uniform values) and approximate matching (either hamming or edit distance)
template <class T, class positionClass> class SearchStrategy;

// Pointer to a partition function
template <class T, class positionClass>
using PartitionPtr = void (SearchStrategy<T, positionClass>::*)(
    const std::string&, std::vector<Substring>&, const int& numParts,
    const int& maxScore, std::vector<SARangePair>& exactMatchRanges) const;

// Pointer to function that starts the index on a particular search
template <class T, class positionClass>
using StartIdxPtr = void (SearchStrategy<T, positionClass>::*)(
    const Search&, const FMOcc<positionClass>&,
    std::vector<FMOcc<positionClass>>&, std::vector<Substring>&,
    const int&) const;

// Pointer to function that starts the index on a particular search for
// strain-free matching
template <class T, class positionClass>
using StartIdxPtrStrainFree = void (SearchStrategy<T, positionClass>::*)(
    const Search&, const FMOcc<positionClass>&,
    std::vector<FMOcc<positionClass>>&, std::vector<Substring>&, const int& idx,
    std::vector<uint32_t>&, std::vector<uint32_t>&) const;

template <class T, class positionClass> class SearchStrategy {
  protected:
    T& index; // pointer to the index of the text that is searched

    // variables for getting info about strategy used
    PartitionStrategy partitionStrategy; // the partitioning strategy
    DistanceMetric distanceMetric;       // which distance metric to use
    std::string name; // the name of the practical search strategy

    // pointers for correct partitioning and correct distance metric
    PartitionPtr<T, positionClass>
        partitionPtr; // pointer to the partition method
    StartIdxPtr<T, positionClass>
        startIdxPtr; // pointer to start method (hamming or
                     // (naive/optimized) edit distance)
    StartIdxPtrStrainFree<T, positionClass>
        startIdxPtrStrainFree; // pointer to start method (naive/optimized edit
                               // distance) for strain-free matching

    // ----------------------------------------------------------------------------
    // CONSTRUCTOR
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param argument, pointer to the bidirectional FM index to use
     * @param p, partition strategy
     * @param edit, true if edit distance should be used, false if hamming
     * distance should be used
     */
    SearchStrategy(T& argument, PartitionStrategy p,
                   DistanceMetric distanceMetric);

    // ----------------------------------------------------------------------------
    // SANITY CHECKS
    // ----------------------------------------------------------------------------

    /**
     * Static function which generates all error patterns with P parts and K
     * errors.
     * @param P the number of parts
     * @param K the number of allowed errors
     * @param patterns vector to store the error patterns in
     */
    static void genErrorPatterns(int P, int K, std::vector<Pattern>& patterns);

    /**
     * Static function to check if a particular search scheme covers all
     * patterns.
     * @param patterns the error patterns to check
     * @param scheme the search scheme to check
     * @param verbose if true the details about which search covers which
     * pattern will be written to stdout
     */
    static bool coversPatterns(const std::vector<Pattern>& patterns,
                               const std::vector<Search>& scheme, bool verbose);

    // ----------------------------------------------------------------------------
    // PARTITIONING
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, either by uniform range or
     * uniform size
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
    void partition(const std::string& pattern, std::vector<Substring>& parts,
                   const int& numParts, const int& maxScore,
                   std::vector<SARangePair>& exactMatchRanges) const;
    /**
     * Calculates the number of parts for a certain max edit distance. This
     * calculation is strategy dependent
     * @param maxED the maximal allowed edit distance for the aligning
     */
    virtual int calculateNumParts(unsigned int maxED) const = 0;

    // Uniform Partitioning

    /**
     * Splits the pattern into numParts parts, such that each part has the
     * same size
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numparts, how many parts are needed
     * @param maxScore, the maximum allowed edit distance
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partitionUniform(const std::string& pattern,
                          std::vector<Substring>& parts, const int& numParts,
                          const int& maxScore,
                          std::vector<SARangePair>& exactMatchRanges) const;

    // Optimal static partitioning

    /**
     * Splits the pattern into numParts parts, such that each search carries
     * the same weight (on average)
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numparts, how many parts are needed
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void
    partitionOptimalStatic(const std::string& pattern,
                           std::vector<Substring>& parts, const int& numParts,
                           const int& maxScore,
                           std::vector<SARangePair>& exactMatchRanges) const;

    /**
     * Helper function for optimal static partitioning. This function
     * creates the optimal static parts
     * @param pattern, the pattern to partition
     * @param parts, empty vector to which the parts are added
     * @param numParts, how many parts there need to be in the partition
     */
    void setParts(const std::string& pattern, std::vector<Substring>& parts,
                  const int& numParts, const int& maxScore) const;

    /**
     * Function that retrieves the begin positions for optimal static
     * partitioning. If derived class does not implement this function then
     * uniform positions are given.
     * @param numparts, how many parts are needed
     * @returns vector with doubles indicating the position (relative to the
     * length of the pattern) where a seed should be placed.
     */
    virtual const std::vector<double> getBegins(const int& numParts,
                                                const int& maxScore) const {
        std::vector<double> b;
        double u = 1.0 / numParts;
        for (int i = 1; i < numParts; i++) {
            b.push_back(i * u);
        }
        return b;
    }

    // Dynamic Partitioning

    /**
     * Splits the pattern into numParts parts, such that each part has
     * (approximately) the same range. The exactMatchRanges are also
     * calculated.
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numparts, how many parts are needed
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partitionDynamic(const std::string& pattern,
                          std::vector<Substring>& parts, const int& numParts,
                          const int& maxScore,
                          std::vector<SARangePair>& exactMatchRanges) const;

    /**
     * Function that retrieves the seeding positions for dynamic partitioning.
     * If derived class does not implement this function then uniform seeds are
     * given.
     * @param numParts, how many parts are needed
     * @returns vector with doubles indicating the position (relative to the
     * length of the pattern) where a seed should be placed.
     */
    virtual const std::vector<double>
    getSeedingPositions(const int& numParts, const int& maxScore) const {

        double u = 1.0 / (numParts - 1);
        std::vector<double> s;
        for (int i = 1; i < numParts - 1; i++) {
            s.push_back(i * u);
        }
        return s;
    }

    /**
     * Helper function for dynamic partitioning. Seeds the parts.
     * @param pattern, the pattern to partition
     * @param parts, empty vector to which the seeds are added
     * @param numparts, how many parts are needed
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     * @returns the number of characters used by the seeding operation
     */
    int seed(const std::string& pattern, std::vector<Substring>& parts,
             const int& numParts, const int& maxScore,
             std::vector<SARangePair>& exactMatchRanges) const;
    /**
     * Function that retrieves the weights for dynamic partitioning.
     * If derived class does not implement this function then uniform weights
     * are given.
     * @param numparts, how many parts are needed
     * @returns vector with weights
     */
    virtual const std::vector<int> getWeights(const int& numParts,
                                              const int& maxScore) const {
        std::vector<int> w(numParts, 1);
        return w;
    }

    /**
     * Helper function for dynamic partitioning. This function extends the parts
     * so that nothing of the pattern is not allocated to any part. This does
     * not keep track of the ranges over the suffix array, so should only be
     * called if this does not matter (e.g. when the parts that can be extended
     * all correspond to empty ranges)
     * @param pattern, the pattern that is split
     * @param parts, the current parts, they are updated so that all characters
     * of pattern are part of exactly one part
     */
    void extendParts(const std::string& pattern,
                     std::vector<Substring>& parts) const;

    // ----------------------------------------------------------------------------
    // (APPROXIMATE) MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Creates all searches for this specific strategy. This is strategy
     * dependent
     * @param maxED the maximal allowed edit distance for the aligning
     */
    virtual const std::vector<Search>&
    createSearches(unsigned int maxED) const = 0;

    /**
     * Executes the search recursively. If U[0] != 1, then the search will
     * start at pi[0], else the search will start with idx i and U[i]!=0 and
     * U[j]=0 with j < i
     * @param s, the search to follow
     * @param parts, the parts of the pattern
     * @param allMatches, vector to add occurrences to
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts
     */
    void doRecSearch(const Search& s, std::vector<Substring>& parts,
                     std::vector<FMOcc<positionClass>>& allMatches,
                     const std::vector<SARangePair>& exactMatchRanges) const;

    /**
     * Starts the index with hamming distance
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     */
    void startIndexHamming(const Search& s,
                           const FMOcc<positionClass>& startMatch,
                           std::vector<FMOcc<positionClass>>& occ,
                           std::vector<Substring>& parts,
                           const int& idx) const {
        index.recApproxMatchHamming(s, startMatch, occ, parts, idx);
    }

    /**
     * Starts the index with edit distance and optimized alignment for the
     * edit distance metric
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     */
    void startIndexEditOptimized(const Search& s,
                                 const FMOcc<positionClass>& startMatch,
                                 std::vector<FMOcc<positionClass>>& occ,
                                 std::vector<Substring>& parts,
                                 const int& idx) const {

        // create unnecessary left and right nodes
        std::vector<uint32_t> leftNodes, rightNodes;
        index.recApproxMatchEditOptimized(s, startMatch, occ, parts, leftNodes,
                                          rightNodes, idx);
    }

    /**
     * Starts the index with edit distance and optimized alignment for the
     * edit distance metric
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     */
    void startIndexEditOptimizedStrainfree(
        const Search& s, const FMOcc<positionClass>& startMatch,
        std::vector<FMOcc<positionClass>>& occ, std::vector<Substring>& parts,
        const int& idx, std::vector<uint32_t>& leftNodes,
        std::vector<uint32_t>& rightNodes) const {

        index.recApproxMatchEditOptimized(s, startMatch, occ, parts, leftNodes,
                                          rightNodes, idx);
    }

    /**
     * Starts the index with naive edit distance (= redundancy between parts
     * of a search)
     *
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     */
    void startIndexEditNaive(const Search& s,
                             const FMOcc<positionClass>& startMatch,
                             std::vector<FMOcc<positionClass>>& occ,
                             std::vector<Substring>& parts,
                             const int& idx) const {

        // create unnecessary left and right nodes
        std::vector<uint32_t> leftNodes, rightNodes;
        index.recApproxMatchEditNaive(s, startMatch, occ, parts, leftNodes,
                                      rightNodes, idx);
    }

    /**
     * Starts the index with naive edit distance (= redundancy between parts
     * of a search)
     *
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     * @param leftNodes the left nodepath (from origin to leftmost point)
     * @param rightNodes the right node path (from origin to rightmost point)
     */
    void startIndexEditNaiveStrainfree(
        const Search& s, const FMOcc<positionClass>& startMatch,
        std::vector<FMOcc<positionClass>>& occ, std::vector<Substring>& parts,
        const int& idx, std::vector<uint32_t>& leftNodes,
        std::vector<uint32_t>& rightNodes) const {

        index.recApproxMatchEditNaive(s, startMatch, occ, parts, leftNodes,
                                      rightNodes, idx);
    }

  public:
    // ----------------------------------------------------------------------------
    // Destructor
    // ----------------------------------------------------------------------------
    virtual ~SearchStrategy() {
    }
    // ----------------------------------------------------------------------------
    // INFORMATION
    // ----------------------------------------------------------------------------

    /**
     * Retrieves the name of this strategy, derived classes should set a
     * meaningful name
     */
    std::string getName() const {
        return name;
    }
    /**
     * Retrieve the  partitioning strategy in string format
     */
    std::string getPartitioningStrategy() const;

    PartitionStrategy getPartitioningStrategyInt() const;

    /**
     * Retrieve the distance metric in string format
     */
    std::string getDistanceMetric() const;

    /**
     * Retrieves the text of the index (for debugging purposes)
     */
    std::string getText() const {
        return index.getText();
    }

    /**
     * Mathces a pattern approximately using this strategy
     * @param pattern, the pattern to match
     * @param maxED, the maximal allowed edit distance (or  hamming
     * distance)
     */
    virtual std::vector<TextOccurrence> matchApprox(const std::string& pattern,
                                                    length_t maxED) const;

    length_t getNodes() const {
        return index.getNodes();
    }

    length_t getMatrixElements() const {
        return index.getMatrixElements();
    }

    length_t getTotal() const {
        return index.getTotalReported();
    }
};

template <class T, class positionClass>
class SearchStrategyDBG : virtual public SearchStrategy<T, positionClass> {
    friend class StrainFreeMapper;

  protected:
    bool strainFree;

    /**
     * Calculates the number of parts for a certain max edit distance. This
     * calculation is strategy dependent
     * @param maxED the maximal allowed edit distance for the aligning
     */
    virtual int calculateNumParts(unsigned int maxED) const = 0;

    /**
     * Creates all searches for this specific strategy. This is strategy
     * dependent
     * @param maxED the maximal allowed edit distance for the aligning
     */
    virtual const std::vector<Search>&
    createSearches(unsigned int maxED) const = 0;

  public:
    SearchStrategyDBG(T& index, PartitionStrategy p = DYNAMIC,
                      DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric),
          strainFree(index.isStrainFree()) {
    }

    virtual std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>
    matchApproxSFI(const std::string& pattern, length_t maxED) const;
};

// ============================================================================
// CLASS SEARCHSTRATEGY
// ============================================================================

// This is a derived class of SearchStrategy. It creates a custom scheme using
// files provided by the user. It takes a specified folder in which a file
// "name.txt", containing the name of the scheme on the first line, and for each
// supported distance score a subfolder exists. Such a subfolder has as name the
// distance score. Each subfolder must contain at least a file "searches.txt".
// The different searches of the scheme for this distance score should be
// written on separate lines of this file. Each search consists out of three
// arrays, pi, L and U, the arrays are separated by a single space. Each array
// is written between curly braces {} and the different values are separated.
// The pi array must be zero-based.
//
// The different subfolders can also contain files for static and dynamic
// partitioning. The file "static_partitioning.txt" should consist out of one
// line of space-separated percentages (between 0 and 1 - exclusive). These
// percentages point to start positions of the second to the last part of the
// partitioned pattern (relative to the size of the pattern). Hence, if a search
// scheme partitions a pattern in k parts, then k+1 percentages should be
// provided. The file "dynamic_partitioning.txt" should consist out of two
// lines. The first line contains k-1 space-separated percentages. These
// percentages are the seeding positions of the middle parts (the first and last
// part are seeded at the begin and end and thus do not need a percentage). Note
// that this line can be empty if the pattern is partitioned into 2 parts. The
// second line contains k integers, where k is the number of parts. Each integer
// corresponds to the weight given to that part in the dynamic partitioning
// process.
template <class T, class positionClass> class CustomSearchStrategy;

// Pointer to the correct getBegins() function for static partitioning
template <class T, class positionClass>
using GetBeginsPtr = const std::vector<double> (
    CustomSearchStrategy<T, positionClass>::*)(const int& numParts,
                                               const int& maxScore) const;

// Pointer to the correct getSeedingPositions() function for dynamic
// partitioning
template <class T, class positionClass>
using GetSeedingPositionsPtr = const std::vector<double> (
    CustomSearchStrategy<T, positionClass>::*)(const int& numParts,
                                               const int& maxScore) const;
// Pointer to the correct getWeights() function for dynamic partitioning
template <class T, class positionClass>
using GetWeightsPtr = const std::vector<int> (
    CustomSearchStrategy<T, positionClass>::*)(const int& numParts,
                                               const int& maxScore) const;

template <class T, class positionClass>
class CustomSearchStrategy : virtual public SearchStrategy<T, positionClass> {

  protected:
    std::vector<std::vector<Search>> schemePerED = {
        {}, {}, {}, {}, {}}; // the search schemes for each distance score, for
                             // now only scores 1 to 4 are supported
    std::vector<bool> supportsMaxScore = {
        false, false, false, false,
        false}; // if a particular distance score is supported

    // static partitioning
    std::vector<std::vector<double>> staticPositions = {
        {}, {}, {}, {}, {}}; // the static positions per score
    std::vector<GetBeginsPtr<T, positionClass>> beginsPointer = {
        &CustomSearchStrategy::getBeginsDefault,
        &CustomSearchStrategy::getBeginsDefault,
        &CustomSearchStrategy::getBeginsDefault,
        &CustomSearchStrategy::getBeginsDefault,
        &CustomSearchStrategy::getBeginsDefault}; // pointer to the correct
                                                  // getBegins() function,
                                                  // either default or custom

    // dynamic partitioning
    std::vector<std::vector<double>> seedingPositions = {
        {}, {}, {}, {}, {}}; // the seeds for dynamic partitioning per score
    std::vector<std::vector<int>> weights = {
        {}, {}, {}, {}, {}}; // the weights for dynamic partitioning per score

    std::vector<GetSeedingPositionsPtr<T, positionClass>> seedingPointer = {
        &CustomSearchStrategy::getSeedingPositionsDefault,
        &CustomSearchStrategy::getSeedingPositionsDefault,
        &CustomSearchStrategy::getSeedingPositionsDefault,
        &CustomSearchStrategy::getSeedingPositionsDefault,
        &CustomSearchStrategy::
            getSeedingPositionsDefault}; // pointer to the correct
                                         // getSeedingPositions() function,
                                         // either default or custom

    std::vector<GetWeightsPtr<T, positionClass>> weightsPointers = {
        &CustomSearchStrategy::getWeightsDefault,
        &CustomSearchStrategy::getWeightsDefault,
        &CustomSearchStrategy::getWeightsDefault,
        &CustomSearchStrategy::getWeightsDefault,
        &CustomSearchStrategy::getWeightsDefault}; // pointer to the correct
                                                   // getWeigths() function,
                                                   // either default or custom

    /**
     * Retrieves the search scheme from a folder, also checks if the scheme is
     * valid
     * @param pathToFolder the path to the folder containing the search scheme
     * @param verbose if the sanity check should be verbose
     */
    void getSearchSchemeFromFolder(std::string pathToFolder, bool verbose);

    /**
     * If the values provided for dynamic partitioning for the given max score
     * are valid (i.e. strictly increasing and between 0 and 1). Will throw a
     * runtime error if this is not the case
     * @param maxScore, the score to check
     */
    void sanityCheckDynamicPartitioning(const int& maxScore) const;

    /**
     * If the values provided for static partitioning for the given max score
     * are valid (i.e. strictly increasing and between 0 and 1). Will throw a
     * runtime error if this is not the case
     * @param maxScore, the score to check
     */
    void sanityCheckStaticPartitioning(const int& maxScore) const;

    /**
     * Parse the search from a line.
     * @param line the line to parse
     * @returns the parsed line as a search, if the line is not valid a runtime
     * error will be thrown.
     */
    Search makeSearch(const std::string& line) const;

    /**
     * Parses an array from a string.
     * @param vectorString the string to parse
     * @param vector the vector with the parsed array as values
     */
    void getVector(const std::string& vectorString,
                   std::vector<int>& vector) const;

    /**
     * Checks wether the connectivity property is satisfied for all searches and
     * if all error patterns are covered for all supported scores. Will throw a
     * runtime error if one of these is not satisfied
     * @param verbose if the information about which search covers which pattern
     * should be written to standard out
     */
    void sanityCheck(bool verbose) const;

    // static partitioning
    /**
     * Gets the static positions in the default manner.
     * @param numParts, the number of parts of the pattern
     * @param maxScore, the maximal allowed score
     */
    const std::vector<double> getBeginsDefault(const int& numParts,
                                               const int& maxScore) const {
        return SearchStrategy<T, positionClass>::getBegins(numParts, maxScore);
    }

    /**
     * Gets the static positions in the custom manner (i.e. with values provided
     * by user in "static_partitioning.txt").
     * @param numParts, the number of parts of the pattern
     * @param maxScore, the maximal allowed score
     */
    const std::vector<double> getBeginsCustom(const int& numParts,
                                              const int& maxScore) const {
        return staticPositions[maxScore - 1];
    }
    /**
     * Overridden function of the base class. Retrieves the begin positions in
     * the custom way if values were provided in a "static_partitioning.txt"
     * file, otherwise the base class function will be called;
     */
    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*beginsPointer[maxScore - 1])(numParts, maxScore);
    }

    // dynamic partitioning
    /**
     * Gets the seeding positions in the default manner.
     * @param numParts, the number of parts of the pattern
     * @param maxScore, the maximal allowed score
     */
    const std::vector<double>
    getSeedingPositionsDefault(const int& numParts, const int& maxScore) const {
        return SearchStrategy<T, positionClass>::getSeedingPositions(numParts,
                                                                     maxScore);
    }

    /**
     * Gets the seeding positions in the custom manner (i.e. with values
     * provided by user in "dynamic_partitioning.txt").
     * @param numParts, the number of parts of the pattern
     * @param maxScore, the maximal allowed score
     */
    const std::vector<double>
    getSeedingPositionsCustom(const int& numParts, const int& maxScore) const {
        return seedingPositions[maxScore - 1];
    }

    /**
     * Overridden function of the base class. Retrieves the seeding positions in
     * the custom way if values were provided in a "dynamic_partitioning.txt"
     * file, otherwise the base class function will be called;
     */
    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*seedingPointer[maxScore - 1])(numParts, maxScore);
    }

    /**
     * Gets the seeding positions in the default manner.
     * @param numParts, the number of parts of the pattern
     * @param maxScore, the maximal allowed score
     */
    const std::vector<int> getWeightsDefault(const int& numParts,
                                             const int& maxScore) const {
        return SearchStrategy<T, positionClass>::getWeights(numParts, maxScore);
    }

    /**
     * Gets the weights in the custom manner (i.e. with values
     * provided by user in "dynamic_partitioning.txt").
     * @param numParts, the number of parts of the pattern
     * @param maxScore, the maximal allowed score
     */
    const std::vector<int> getWeightsCustom(const int& numParts,
                                            const int& maxScore) const {
        return weights[maxScore - 1];
    }
    /**
     * Overridden function of the base class. Retrieves the weights positions in
     * the custom way if values were provided in a "dynamic_partitioning.txt"
     * file, otherwise the base class function will be called;
     */
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*weightsPointers[maxScore - 1])(numParts, maxScore);
    }

  public:
    CustomSearchStrategy(T& index, std::string pathToFolder,
                         PartitionStrategy p = DYNAMIC,
                         DistanceMetric metric = EDITOPTIMIZED,
                         bool verbose = false)
        : SearchStrategy<T, positionClass>(index, p, metric) {

        getSearchSchemeFromFolder(pathToFolder, verbose);
    }

    int calculateNumParts(unsigned int maxED) const {
        assert(supportsMaxScore[maxED - 1]);
        return schemePerED[maxED - 1][0].getNumParts();
    }
    const std::vector<Search>& createSearches(unsigned int maxED) const {
        assert(supportsMaxScore[maxED - 1]);
        return schemePerED[maxED - 1];
    }
};

template <class T, class positionClass>
class CustomSearchStrategyDBG : public SearchStrategyDBG<T, positionClass>,
                                public CustomSearchStrategy<T, positionClass> {
  private:
    int calculateNumParts(unsigned int maxED) const {
        return CustomSearchStrategy<T, positionClass>::calculateNumParts(maxED);
    }

    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return CustomSearchStrategy<T, positionClass>::createSearches(maxED);
    }

  public:
    CustomSearchStrategyDBG(T& index, std::string pathToFolder,
                            PartitionStrategy p = DYNAMIC,
                            DistanceMetric metric = EDITOPTIMIZED,
                            bool verbose = false)
        : SearchStrategy<T, positionClass>(index, p, metric),
          SearchStrategyDBG<T, positionClass>(index, p, metric),
          CustomSearchStrategy<T, positionClass>(index, pathToFolder, p, metric,
                                                 verbose) {
    }
};

// ============================================================================
// CLASS NaiveBackTrackingStrategy
// ============================================================================

// Matches a pattern using the naive backtracking strategy.
template <class T, class positionClass>
class NaiveBackTrackingStrategy
    : virtual public SearchStrategy<T, positionClass> {
  protected:
    std::vector<Search> searches = {};
    int calculateNumParts(unsigned int maxED) const {
        return 1;
    }
    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return searches;
    }

  public:
    virtual std::vector<TextOccurrence> matchApprox(const std::string& pattern,
                                                    length_t maxED) const {

        SearchStrategy<T, positionClass>::index.resetCounters();
        if (maxED == 0) {

            auto result =
                SearchStrategy<T, positionClass>::index.exactMatches(pattern);
            std::vector<TextOccurrence> returnvalue;
            for (length_t startpos : result) {
                returnvalue.emplace_back(
                    Range(startpos, startpos + pattern.size()), 0);
                returnvalue.back().generateOutput();
            }
            return returnvalue;
        }
        return SearchStrategy<T, positionClass>::index.approxMatchesNaive(
            pattern, maxED);
    }

    NaiveBackTrackingStrategy(T& index, PartitionStrategy p = DYNAMIC,
                              DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric) {
        SearchStrategy<T, positionClass>::name = "Naive backtracking";
    };
};

template <class T, class positionClass>
class NaiveBackTrackingStrategyDBG
    : public SearchStrategyDBG<T, positionClass>,
      public NaiveBackTrackingStrategy<T, positionClass> {
  private:
    int calculateNumParts(unsigned int maxED) const {
        return NaiveBackTrackingStrategy<T, positionClass>::calculateNumParts(
            maxED);
    }

    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return NaiveBackTrackingStrategy<T, positionClass>::createSearches(
            maxED);
    }

  public:
    virtual std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>
    matchApproxSFI(const std::string& pattern, length_t maxED) const {

        SearchStrategy<T, positionClass>::index.resetCounters();
        if (maxED == 0) {

            auto result =
                SearchStrategy<T, positionClass>::index.ExactMatchSFI(pattern);
            std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>
                paths;
            for (auto occ : result) {
                occ.generateOutput();
                paths[occ.getNodePath()].emplace_back(occ);
            }
            return paths;
        }
        return SearchStrategy<T, positionClass>::index.approxMatchesNaiveSFI(
            pattern, maxED);
    }

    NaiveBackTrackingStrategyDBG(T& index, PartitionStrategy p = DYNAMIC,
                                 DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric),
          SearchStrategyDBG<T, positionClass>(index, p, metric),
          NaiveBackTrackingStrategy<T, positionClass>(index, p, metric) {
    }
};

// ============================================================================
// HARDCODED CUSTOM CLASSES
// ============================================================================

template <class T, class positionClass>
class KucherovKplus1 : virtual public SearchStrategy<T, positionClass> {
  protected:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 0}, {0, 1}),
        Search::makeSearch({1, 0}, {0, 0}, {0, 1})};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}),
        Search::makeSearch({2, 1, 0}, {0, 0, 0}, {0, 1, 2}),
        Search::makeSearch({1, 0, 2}, {0, 0, 1}, {0, 1, 2})};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 3, 3}),
        Search::makeSearch({1, 0, 2, 3}, {0, 0, 1, 1}, {0, 1, 3, 3}),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 1, 3, 3}),
        Search::makeSearch({3, 2, 1, 0}, {0, 0, 1, 1}, {0, 1, 3, 3})};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 2, 2, 4, 4}),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 4, 4}),
        Search::makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 4, 4}),
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 4, 4}),
        Search::makeSearch({3, 2, 4, 1, 0}, {0, 0, 0, 1, 1}, {0, 1, 2, 4, 4}),
        Search::makeSearch({2, 1, 0, 3, 4}, {0, 0, 0, 1, 3}, {0, 1, 2, 4, 4}),
        Search::makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 2, 4}, {0, 1, 2, 4, 4}),
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 3, 4}, {0, 0, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    const std::vector<std::vector<double>> seedingPositions = {
        {}, {0.57}, {0.38, 0.65}, {0.38, 0.55, 0.73}};

    const std::vector<std::vector<int>> weights = {
        {1, 1}, {39, 10, 40}, {400, 4, 5, 400}, {100, 5, 1, 6, 105}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.41, 0.7}, {0.25, 0.50, 0.75}, {0.27, 0.47, 0.62, 0.81}};
    int calculateNumParts(unsigned int maxED) const {
        return maxED + 1;
    }
    const std::vector<Search>& createSearches(unsigned int maxED) const {
        assert(maxED >= 1);
        assert(maxED <= 4);

        return schemePerED[maxED - 1];
    }
    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

  public:
    KucherovKplus1(T& index, PartitionStrategy p = DYNAMIC,
                   DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric) {
        SearchStrategy<T, positionClass>::name = "KUCHEROV K + 1";
    };
};

template <class T, class positionClass>
class KucherovKplus1DBG : public SearchStrategyDBG<T, positionClass>,
                          public KucherovKplus1<T, positionClass> {
  private:
    int calculateNumParts(unsigned int maxED) const {
        return KucherovKplus1<T, positionClass>::calculateNumParts(maxED);
    }

    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return KucherovKplus1<T, positionClass>::createSearches(maxED);
    }

  public:
    KucherovKplus1DBG(T& index, PartitionStrategy p = DYNAMIC,
                      DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric),
          SearchStrategyDBG<T, positionClass>(index, p, metric),
          KucherovKplus1<T, positionClass>(index, p, metric) {
    }
};

template <class T, class positionClass>
class KucherovKplus2 : virtual public SearchStrategy<T, positionClass> {
  protected:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}),
        Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1})};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 1, 2}),
        Search::makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        Search::makeSearch({1, 2, 3, 0}, {0, 0, 0, 1}, {0, 0, 1, 2}),
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 2}, {0, 0, 2, 2})};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 2, 3, 3}),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 2, 2, 3}),
        Search::makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 1}, {0, 1, 1, 3, 3}),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 1, 2}, {0, 0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 4}),
        Search::makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 4}),
        Search::makeSearch({5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 2, 4, 4}),
        Search::makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 1, 2},
                           {0, 1, 1, 3, 4, 4}),
        Search::makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 2, 3},
                           {0, 1, 1, 2, 4, 4}),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 1, 3, 3},
                           {0, 0, 3, 3, 4, 4}),
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 3, 3, 3},
                           {0, 0, 3, 3, 4, 4}),
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 4, 4},
                           {0, 0, 2, 4, 4, 4}),
        Search::makeSearch({2, 3, 1, 0, 4, 5}, {0, 0, 0, 1, 2, 4},
                           {0, 0, 2, 2, 4, 4}),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 4, 4},
                           {0, 0, 1, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    int calculateNumParts(unsigned int maxED) const {
        return maxED + 2;
    }
    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.48, 0.55}, {0.4, 0.63, 0.9}, {0.34, 0.5, 0.65, 0.7}};

    const std::vector<std::vector<int>> weights = {{11, 10, 1},
                                                   {400, 4, 1, 800},
                                                   {6, 3, 2, 1, 1},
                                                   {52, 42, 16, 14, 1, 800}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.47, 0.94},
        {0.35, 0.50, 0.65},
        {0.22, 0.44, 0.66, 0.88},
        {0.18, 0.37, 0.53, 0.69, 0.83}};

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

  public:
    KucherovKplus2(T& index, PartitionStrategy p = DYNAMIC,
                   DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric) {
        SearchStrategy<T, positionClass>::name = "KUCHEROV K + 2";
    };
};

template <class T, class positionClass>
class KucherovKplus2DBG : public SearchStrategyDBG<T, positionClass>,
                          public KucherovKplus2<T, positionClass> {
  private:
    int calculateNumParts(unsigned int maxED) const {
        return KucherovKplus2<T, positionClass>::calculateNumParts(maxED);
    }

    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return KucherovKplus2<T, positionClass>::createSearches(maxED);
    }

  public:
    KucherovKplus2DBG(T& index, PartitionStrategy p = DYNAMIC,
                      DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric),
          SearchStrategyDBG<T, positionClass>(index, p, metric),
          KucherovKplus2<T, positionClass>(index, p, metric) {
    }
};

template <class T, class positionClass>
class OptimalKianfar : virtual public SearchStrategy<T, positionClass> {
  protected:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 0}, {0, 1}),
        Search::makeSearch({1, 0}, {0, 1}, {0, 1})};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 2}, {0, 1, 2}),
        Search::makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2}),
        Search::makeSearch({1, 2, 0}, {0, 1, 1}, {0, 1, 2})};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 3}, {0, 2, 3, 3}),
        Search::makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {1, 2, 3, 3}),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 2, 2}, {0, 0, 3, 3})};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 4}, {0, 3, 3, 4, 4}),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {2, 2, 3, 3, 4}),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 3, 3}, {0, 0, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    const std::vector<std::vector<double>> seedingPositions = {
        {}, {0.50}, {0.34, 0.66}, {0.42, 0.56, 0.67}};

    const std::vector<std::vector<int>> weights = {
        {1, 1}, {10, 1, 5}, {1, 1, 1, 1}, {7, 2, 1, 3, 5}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.30, 0.60}, {0.17, 0.69, 0.96}, {0.2, 0.5, 0.6, 0.8}};
    int calculateNumParts(unsigned int maxED) const {
        return maxED + 1;
    }
    const std::vector<Search>& createSearches(unsigned int maxED) const {
        if (maxED < 1 || maxED > 5) {
            throw std::invalid_argument("max ED should be between 1 and 4");
        }
        return schemePerED[maxED - 1];
    }

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

  public:
    OptimalKianfar(T& index, PartitionStrategy p = DYNAMIC,
                   DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric) {
        SearchStrategy<T, positionClass>::name = "OPTIMAL KIANFAR";
    };
};

template <class T, class positionClass>
class OptimalKianfarDBG : public SearchStrategyDBG<T, positionClass>,
                          public OptimalKianfar<T, positionClass> {
  private:
    int calculateNumParts(unsigned int maxED) const {
        return OptimalKianfar<T, positionClass>::calculateNumParts(maxED);
    }

    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return OptimalKianfar<T, positionClass>::createSearches(maxED);
    }

  public:
    OptimalKianfarDBG(T& index, PartitionStrategy p = DYNAMIC,
                      DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric),
          SearchStrategyDBG<T, positionClass>(index, p, metric),
          OptimalKianfar<T, positionClass>(index, p, metric) {
    }
};

// ============================================================================
// CLASS O1StarSearchStrategy
// ============================================================================

// A concrete derived class of SearchStrategy. The strategy here is founded
// on this observation: if x errors are allowed and the pattern is divided
// up in (x
// + 2) parts then every match with max x erros contains a seed consisting
// of n parts, where the first and last part of the seed contain no errors
// and all parts inbetween these contain exacly one error. (2 <= n <= x +
// 2)
template <class T, class positionClass>
class O1StarSearchStrategy : virtual public SearchStrategy<T, positionClass> {
  protected:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}),
        Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1})};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        Search::makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 0, 2, 2})};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        Search::makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}),
        Search::makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}),
        Search::makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}),
        Search::makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 0, 4, 4, 4, 4}),
    };

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    int calculateNumParts(unsigned int maxED) const {
        return maxED + 2;
    }
    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.51, 0.93}, {0.34, 0.64, 0.88}, {0.28, 0.48, 0.63, 0.94}};

    const std::vector<std::vector<int>> weights = {
        {11, 10, 1}, {20, 11, 11, 10}, {3, 2, 2, 1, 1}, {1, 2, 2, 1, 2, 1}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.50, 0.96},
        {0.26, 0.64, 0.83},
        {0.22, 0.46, 0.67, 0.95},
        {0.19, 0.37, 0.57, 0.74, 0.96}};

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

  public:
    O1StarSearchStrategy(T& index, PartitionStrategy p = DYNAMIC,
                         DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric) {
        SearchStrategy<T, positionClass>::name = "01*0";
    };
};

template <class T, class positionClass>
class O1StarSearchStrategyDBG : public SearchStrategyDBG<T, positionClass>,
                                public O1StarSearchStrategy<T, positionClass> {
  private:
    int calculateNumParts(unsigned int maxED) const {
        return O1StarSearchStrategy<T, positionClass>::calculateNumParts(maxED);
    }

    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return O1StarSearchStrategy<T, positionClass>::createSearches(maxED);
    }

  public:
    O1StarSearchStrategyDBG(T& index, PartitionStrategy p = DYNAMIC,
                            DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric),
          SearchStrategyDBG<T, positionClass>(index, p, metric),
          O1StarSearchStrategy<T, positionClass>(index, p, metric) {
    }
};

template <class T, class positionClass>
class ManBestStrategy : virtual public SearchStrategy<T, positionClass> {
  protected:
    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 4},
                           {0, 3, 3, 3, 4, 4}),
        Search::makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 2, 2, 3, 3, 4}),
        Search::makeSearch({2, 1, 3, 4, 5, 0}, {0, 1, 1, 1, 1, 1},
                           {0, 2, 2, 3, 3, 4}),
        Search::makeSearch({3, 2, 1, 4, 5, 0}, {0, 1, 2, 2, 2, 2},
                           {0, 1, 2, 3, 3, 4}),
        Search::makeSearch({5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 3, 3},
                           {0, 0, 4, 4, 4, 4})};

    int calculateNumParts(unsigned int maxED) const {
        return maxED + 2;
    }
    const std::vector<Search>& createSearches(unsigned int maxED) const {
        assert(maxED == 4);
        return ED4;
    }

    const std::vector<double> seedingPositions = {0.35, 0.59, 0.67, 0.9};

    const std::vector<int> weights = {89, 15, 90, 1, 48, 84};

    const std::vector<double> staticPositions = {0.24, 0.43, 0.62, 0.73, 0.77};

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        assert(maxScore == 4);
        return staticPositions;
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        assert(maxScore == 4);
        return weights;
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        assert(maxScore == 4);
        return seedingPositions;
    }

  public:
    ManBestStrategy(T& index, PartitionStrategy p = DYNAMIC,
                    DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric) {
        SearchStrategy<T, positionClass>::name = "MANBEST";
    };
};

template <class T, class positionClass>
class ManBestStrategyDBG : public SearchStrategyDBG<T, positionClass>,
                           public ManBestStrategy<T, positionClass> {
  private:
    int calculateNumParts(unsigned int maxED) const {
        return ManBestStrategy<T, positionClass>::calculateNumParts(maxED);
    }

    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return ManBestStrategy<T, positionClass>::createSearches(maxED);
    }

  public:
    ManBestStrategyDBG(T& index, PartitionStrategy p = DYNAMIC,
                       DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric),
          SearchStrategyDBG<T, positionClass>(index, p, metric),
          ManBestStrategy<T, positionClass>(index, p, metric) {
    }
};
// ============================================================================
// CLASS PIGEONHOLESEARCHSTRATEGY
// ============================================================================

// A concrete derived class of SearchStrategy. The strategy here is founded
// on this observation: if x errors are allowed and the pattern is divided
// up in (x
// + 1) sections then every approximate match has an exact match with at
// least one of the sections. The strategy iterates over the sections, it
// tries to exactly match the current section, then approximately match the
// pattern before this section and after the pattern after this section
// with the remaining edit distance.
template <class T, class positionClass>
class PigeonHoleSearchStrategy
    : virtual public SearchStrategy<T, positionClass> {

  protected:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 0}, {0, 1}),
        Search::makeSearch({1, 0}, {0, 0}, {0, 1})};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}),
        Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 2, 2}),
        Search::makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2})};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}),
        Search::makeSearch({1, 0, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3}),
        Search::makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        Search::makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    int calculateNumParts(unsigned int maxED) const {
        return maxED + 1;
    }
    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return schemePerED[maxED - 1];
    }

  public:
    PigeonHoleSearchStrategy(T& index, PartitionStrategy p = DYNAMIC,
                             DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric) {
        SearchStrategy<T, positionClass>::name = "PIGEONHOLE";
    };
};

template <class T, class positionClass>
class PigeonHoleSearchStrategyDBG
    : public SearchStrategyDBG<T, positionClass>,
      public PigeonHoleSearchStrategy<T, positionClass> {
  private:
    int calculateNumParts(unsigned int maxED) const {
        return PigeonHoleSearchStrategy<T, positionClass>::calculateNumParts(
            maxED);
    }

    const std::vector<Search>& createSearches(unsigned int maxED) const {
        return PigeonHoleSearchStrategy<T, positionClass>::createSearches(
            maxED);
    }

  public:
    PigeonHoleSearchStrategyDBG(T& index, PartitionStrategy p = DYNAMIC,
                                DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy<T, positionClass>(index, p, metric),
          SearchStrategyDBG<T, positionClass>(index, p, metric),
          PigeonHoleSearchStrategy<T, positionClass>(index, p, metric) {
    }
};

#endif
