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

#include "fmindex.h"
#include "globalParameters.h"
#include "mappingpair.h"
#include "node.h"
#include "rankinterface.h"
#include "rankselectinterface.h"

#include <map>
#include <queue>

// ============================================================================
// CLASS FMINDEXDBG
// ============================================================================

template <class positionClass>
class FMIndexDBG : public FMIndex<positionClass> {
    friend class FMPos;
    friend class FMPosSFR;

  private:
    // De Bruijn parameter
    uint k;
    // Number of strains in the pan-gemome
    uint numberOfStrains;
    // Number of nodes in the graph;
    uint numberOfGraphNodes;

    // Bit vectors for the implicit representation of the compressed de Bruijn
    // graph

    // Set the bits to one for the last entry in the range in the SA for the
    // rightmost k-mer of every node. This bit vector supports rank and select.
    RankSelectInterface B_right;
    // Set the bits to one for the last entry in the range in the reverse SA for
    // the leftmost k-mer of every node. This bit vector supports rank.
    RankInterface B_left;
    // Set the bits to one for every entry in the range in the SA for the
    // rightmost k-mer of every node. This bit vector supports rank.
    RankInterface B_right_full;

    // Mapping of right node IDs
    std::vector<MappingPair> mapping_right;
    // Mapping of left node IDs
    std::vector<int> mapping_left;
    // List of sorted startpositions of the strains along with the position of
    // the sentinel character
    std::vector<length_t> sorted_startpositions{0};
    // Compacted de Bruijn graph: vector of nodes
    std::vector<Node> G;

    // Counter for the number of visited graph nodes
    thread_local static length_t nodeDBGCounter;
    // Counter for the number of special cases in the linear strain-free
    // filtering process
    thread_local static length_t filterSpecialCaseCounter;

    // the (sparse) reverse suffix array of the reference genome, this is only
    // initialized in the build process
    SparseSuffixArray sparseRevSA;

    // Boolean that stores the filtering option in case of strain-free matching:
    // linear or complete
    bool filteringOptionComplete = false;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR THE BUILDING PROCESS
    // ----------------------------------------------------------------------------

    /**
     * @brief Build a new implicit de Bruijn graph, along with its underlying
     * bidirectional FM-index.
     *
     * @param baseFile Base filename for FM-index
     * @param k k-mer size
     * @param sa_sparse Suffix array sparseness factor
     * @param checkpoint_sparseness Sparseness factors for the checkpoints that
     * aid in finding node identifiers.
     * @param progress Prints progress if true
     * @param option Select algorithm option
     */
    FMIndexDBG(const std::string& baseFile, const int k,
               const std::vector<int>& sa_sparse,
               const std::vector<int>& checkpoint_sparseness,
               const bool progress,
               const SelectOption& option = SelectOption::SIMPLE)
        : FMIndex<positionClass>(0, baseFile), k(k) {

        createFMIndex(baseFile, sa_sparse);

        k_DBG = k;

        this->numberOfSeparationCharacters = 2;

        // Count the number of strains
        numberOfStrains =
            (std::count(this->bwt.begin(), this->bwt.end(), '%') + 1);

        this->numberOfSeparationCharacters = 2;
        // Initialize the bit vectors
        B_right.setN(this->textLength + 1);
        B_right.setOption(RankSelectOption::RANK9SELECT);
        B_left.setN(this->textLength + 1);
        B_left.setOption(RankOption::RANK9);
        B_right_full.setN(this->textLength + 1);
        B_right_full.setOption(RankOption::RANK9);

        std::vector<RankSelectInterface> B_rights(checkpoint_sparseness.size());
        std::vector<RankInterface> B_right_fulls(checkpoint_sparseness.size());

        for (uint j = 0; j < checkpoint_sparseness.size(); j++) {
            B_rights[j].setN(this->textLength + 1);
            B_rights[j].setOption(RankSelectOption::RANK9SELECT);
            B_right_fulls[j].setN(this->textLength + 1);
            B_right_fulls[j].setOption(RankOption::RANK9);
        }

        std::vector<std::vector<MappingPair>> mapping_rights(
            checkpoint_sparseness.size());

        // Create the implicit compressed De Bruijn graph
        buildCompressedGraph(checkpoint_sparseness, progress, B_rights,
                             B_right_fulls, mapping_rights);

        // Write the bitvectors to output files
        {
            for (uint j = 0; j < checkpoint_sparseness.size(); j++) {
                std::string suffix =
                    checkpoint_sparseness[j] == INT32_MAX
                        ? "none"
                        : std::to_string(checkpoint_sparseness[j]);
                std::ofstream ofs(baseFile + ".B.right." + suffix);
                if (!ofs) {
                    throw std::runtime_error("Cannot open file: " + baseFile +
                                             ".B.right." + suffix);
                }
                B_rights[j].write(ofs);
                ofs.close();
            }
        }
        {
            std::ofstream ofs(baseFile + ".B.left");
            if (!ofs) {
                throw std::runtime_error("Cannot open file: " + baseFile +
                                         ".B.left");
            }
            B_left.write(ofs);
            ofs.close();
        }
        // Write the bitvectors to output files
        {
            for (size_t j = 0; j < checkpoint_sparseness.size(); j++) {
                std::string suffix =
                    checkpoint_sparseness[j] == INT32_MAX
                        ? "none"
                        : std::to_string(checkpoint_sparseness[j]);
                std::ofstream ofs(baseFile + ".B.right.full." + suffix);
                if (!ofs) {
                    throw std::runtime_error("Cannot open file: " + baseFile +
                                             ".B.right.full." + suffix);
                }
                B_right_fulls[j].write(ofs);
                ofs.close();
            }
        }

        // Write the compressed De Bruijn graph to an output file
        {
            std::ofstream ofs(baseFile + ".DBG");
            ofs.write((char*)&k, sizeof(k));
            for (std::vector<Node>::iterator itr = G.begin(); itr != G.end();
                 ++itr) {
                itr->write(ofs);
            }
            ofs.close();
        }

        // Write the node ID mappings to an output file
        {
            for (uint j = 0; j < checkpoint_sparseness.size(); j++) {
                std::string suffix =
                    checkpoint_sparseness[j] == INT32_MAX
                        ? "none"
                        : std::to_string(checkpoint_sparseness[j]);
                std::ofstream ofs(baseFile + ".right.map." + suffix,
                                  std::ios::binary);
                if (!ofs) {
                    throw std::runtime_error("Cannot open file: " + baseFile +
                                             ".right.map." + suffix);
                }
                for (std::vector<MappingPair>::iterator itr =
                         mapping_rights[j].begin();
                     itr != mapping_rights[j].end(); ++itr) {
                    itr->write(ofs);
                }
                ofs.close();
            }
        }
        {
            std::ofstream ofs(baseFile + ".left.map", std::ios::binary);
            if (!ofs) {
                throw std::runtime_error("Cannot open file: " + baseFile +
                                         ".left.map");
            }
            ofs.write((char*)&mapping_left[0],
                      mapping_left.size() * sizeof(int));
            ofs.close();
        }
    }

    /**
     * @brief Build the underlying bidirectional FM-index that will be the base
     * of the implicit representation of the compressed de Bruijn graph. This
     * function also writes out the necessary data structures to files for later
     * use.
     *
     * @param baseFN Base filename for FM-index
     * @param sparse_sa Suffix array sparseness factor
     */
    void createFMIndex(const std::string& baseFN,
                       const std::vector<int>& sparse_sa);

    /**
     * @brief Find the entry in the reverse suffix array of this index. This is
     * computed from the sparse reverse suffix array and the reverse bwt. This
     * function only works when the reverse SA is stored, which is only during
     * the build process.
     *
     * @param index the index to find the entry in the reverse SA of
     * @return length_t - the entry in the SA of the index
     */
    length_t findRevSA(length_t index) const;

    /**
     * @brief Build the compacted longest common prefix array. An entry of 0
     * means that the LCP is shorter than k, an entry of 1 means that the LCP
     * has length k and an entry of 2 means that the LCP is longer than k.
     *
     * @param LCP A 2-bitvector object in which the compacted LCP will be
     * stored. It will be overwritten.
     * @param progress Prints progress if true
     */
    void computeLCP(Bitvec2& LCP, bool progress);

    /**
     * @brief Compute the Br_right and Bl_right bitvectors
     *
     * @param Q Node queue
     */

    /**
     * @brief Build the Br_right and Bl_right bitvectors
     *
     * @param Q the node queue
     * @param Br_right bit vector in which Br_right will be stored, will be
     * overwritten
     * @param Bl_right bit vector in which Br_right will be stored, will be
     * overwritten
     * @param progress Prints progress if true
     */
    void computeBitVectors(std::queue<int>& Q, Bitvec& Br_right,
                           Bitvec& Bl_right, bool progress);

    /**
     * @brief For an ??-interval [i,j[, the function call getIntervals(i,j)
     * returns the vector of all c??-intervals.
     *
     * @param i Left bound of ??-interval
     * @param j Right bound of ??-interval
     * @return std::vector<std::pair<char, Range>> - vector of all c??-intervals
     */
    std::vector<std::pair<char, Range>> getIntervals(length_t i, length_t j);

    /**
     * @brief Find the range over the reverse suffix array that corresponds to
     * the first k-mer of the suffix at the indexInSA'th position in the SA
     *
     * @param indexInSA index in the suffix array for the suffix of which the
     * first k-mer is of interest
     * @param isEndNode indicates whether the node we are looking at is an end
     * node, since these have a different convention
     * @return Range - the range over the reverse suffix array corresponding to
     * the k-mer of interest
     */
    Range getReverseRange(length_t indexInSA, bool isEndNode);

    /**
     * @brief Set the edge mapping of a node. The edge mapping maps the ranks of
     * the edges in the SA to their ranks in the reverse SA.
     *
     * @param id The identifier of the node for which the mapping must be
     * constructed.
     */
    void setEdgeMapping(length_t id);

    /**
     * @brief Build the implicit compressed de Bruijn graph
     *
     * @param checkpoint_sparseness Sparseness factors for the checkpoints that
     * aid in finding node identifiers.
     * @param progress Prints progress if true
     * @param B_rights B_right bit vectors for different checkpoint sparseness
     * factors
     * @param B_right_fulls B_right_full bit vectors for different checkpoint
     * sparseness factors
     * @param mapping_rights mapping_right mappings for different checkpoint
     * sparseness factors
     */
    void
    buildCompressedGraph(const std::vector<int>& checkpoint_sparseness,
                         bool progress,
                         std::vector<RankSelectInterface>& B_rights,
                         std::vector<RankInterface>& B_right_fulls,
                         std::vector<std::vector<MappingPair>>& mapping_rights);

    // ----------------------------------------------------------------------------
    // HELPER ROUTINES FOR MAPPING
    // ----------------------------------------------------------------------------

    /**
     * @brief If i is an index for the suffix array, find the node in the
     * compressed de Bruijn graph that corresponds to T[SA[i]..SA[i]+k[. The
     * function assumes that the given k-mer is the last k-mer of the node.
     *
     * @param i Index in the suffix array indicating the last k-mer of a node
     * @return int - the ID of the corresponding node in the compressed de
     * Bruijn graph
     */
    int findIDLast(length_t i) const;

    /**
     * @brief If i is an index for the suffix array, find the node u in the
     * compressed de Bruijn graph that corresponds to T[SA[i]..SA[i]+k[. The
     * function assumes that the given k-mer is the first k-mer of the node.
     *
     * @param i Index in the suffix array indicating the first k-mer of a node
     * @return int - the ID of the corresponding node in the compressed de
     * Bruijn graph
     */
    int findIDFirst(length_t i) const;

    /**
     * @brief If i is an index for the suffix array, find the node in the
     * compressed de Bruijn graph that corresponds to T[SA[i]..SA[i]+k[.
     * T[SA[j-1]..SA[j-1]+k[ can be at any location in the node.
     *
     * @param i Index in the suffix array indicating a certain k-mer of a
     * node
     * @param id the ID of the corresponding node in the compressed de Bruijn
     * graph (to be filled in)
     * @param l the number of characters that are before the k-mer in the node
     * (to be filled in)
     */

    /**
     * @brief If i is an index for the suffix array, find the node in the
     * compressed de Bruijn graph that corresponds to T[SA[i]..SA[i]+k[.
     * T[SA[j-1]..SA[j-1]+k[ can be at any location in the node.
     *
     * @param i Index in the suffix array indicating a certain k-mer of a
     * node
     * @param id the ID of the corresponding node in the compressed de Bruijn
     * graph (to be filled in)
     * @param l the number of characters that are before the k-mer in the node
     * (to be filled in)
     * @param offset the offset of the edge in the reverse SA corresponding to
     * index i in the regular SA (to be filled in)
     */
    void findIDandOffset(length_t i, uint32_t& id, uint32_t& l,
                         uint32_t& offset) const;

    /**
     * @brief Find the identifier of the successor of node id by following the
     * 'offset'th edge (with respect to the reverse SA) through node id. Also
     * find the new offset (with respect to the reverse SA) in the new node.
     *
     * @param id Identifier of the current node
     * @param offset_reverse Offset of the edge to follow to the next node (with
     * respect to the reverse SA). Will be updated.
     * @return int - The identifier of the successor of interest
     */
    int jumpToSuccessorThroughEdge(uint32_t id, uint32_t& offset_reverse) const;

    /**
     * @brief Find the identifier of the predecessor of node id by following the
     * 'offset'th edge (with respect to the regular SA) through node id. Also
     * find the new offset (with respect to the regular SA) in the new node.
     *
     * @param id Identifier of the current node
     * @param offset Offset of the edge to follow to the previous node (with
     * respect to the regular SA). Will be updated.
     * @return int - The identifier of the predecessor of interest
     */
    int jumpToPredecessorThroughEdge(uint32_t id, uint32_t& offset) const;

    /**
     * @brief Find the identifier of the successor of node id by appending a
     * character to its substring.
     *
     * @param id Identifier of the current node
     * @param id_successor Identifier of the successor (to be filled in)
     * @param posInAlphabet Position in the alfabet of the character to add
     * @param reverse_offset The offset in the new interval (with respect to the
     * reverse SA). This is needed for end nodes.
     * @return true - a successor was found
     * @return false - no successor was found
     */
    bool jumpToSuccessorWithChar(uint32_t id, uint32_t& id_successor,
                                 uint32_t posInAlphabet,
                                 uint32_t reverse_offset = 0) const;

    /**
     * @brief Find the identifier of the predecessor of node id by prepending a
     * character to its substring.
     *
     * @param id Identifier of the current node
     * @param id_predecessor Identifier of the predecessor (to be filled in)
     * @param posInAlphabet Position in the alfabet of the character to add
     * @param offset Placeholder variable such that jumpToPredecessorWithChar
     * and jumpToSuccessorWithChar have the same parameters. This variable is
     * not used.
     * @return true - a predecessor was found
     * @return false - no predecessor was found
     */
    bool jumpToPredecessorWithChar(uint32_t id, uint32_t& id_predecessor,
                                   uint32_t posInAlphabet,
                                   uint32_t offset = 0) const;

    /**
     * @brief This function maps matches in the SA along with their node path to
     * matches in the text. It takes the ranges of the matches and together with
     * the depth this is matched to a range in the text (this new range has a
     * width of depth). The edit distance is also mapped to this new range in
     * the text. This function also filters out the redundant matches using the
     * maximal allowed edit distance.
     *
     * @param saMatch The approximate match with its ranges in the SA and revSA
     * and its node path that needs to be mapped to a match in the text
     * @return std::vector<TextOccurrenceSFI> - a vector of matches in the text
     * containing a range and the edit distance
     */
    std::vector<TextOccurrenceSFI>
    convertToMatchesInTextSFI(const FMOccSFI<positionClass>& saMatch);

    /**
     * @brief This function maps matches in the SA along with their node path to
     * matches in the text. It takes the ranges of the matches and together with
     * the depth this is matched to a range in the text (this new range has a
     * width of depth). The edit distance is also mapped to this new range in
     * the text. This function also filters out the redundant matches using the
     * maximal allowed edit distance
     *
     * @param ranges the ranges in the SA and revSA of the approximate match
     * @param nodepath the node path of the approximate match
     * @param patternLength the length of the approximate match
     * @param distance the error distance of the approximate match
     * @param shift the shit of the approximate match
     * @return std::vector<TextOccurrenceSFI> - a vector of matches in the text
     * containing a range and the edit distance
     */
    std::vector<TextOccurrenceSFI>
    convertToMatchesInTextSFI(const SARangePair& ranges,
                              const std::vector<uint32_t>& nodepath,
                              const int& patternLength, const int& distance = 0,
                              const length_t& shift = 0);

    /**
     * @brief Check if one of the children of a position is a separation
     * character.
     *
     * @param pos The position of which the children need to be checked
     * @return true if one of the children of the position is a separation
     * character
     * @return false otherwise
     */
    const bool separationIsNext(positionClass pos) const;

    /**
     * @brief Pushes all the children corresponding to the node with ranges
     * equal to parentRanges onto the stack
     *
     * @param parentRanges the ranges to get the children of
     * @param stack the stack to push the children on
     * @param row the row of the parentNode
     * @param trueDepth the true depth of the total match corresponding to
     * parentNode. Defaults to -1
     */
    void extendFMPos(const SARangePair& parentRanges,
                     std::vector<FMPosExt<positionClass>>& stack, int row,
                     int trueDepth = -1) override;

    /**
     * @brief Pushes all the children corresponding to the this position onto
     * the stack
     *
     * @param pos the position to get the children of
     * @param stack the stack to push the children on
     */
    void extendFMPos(const positionClass& pos,
                     std::vector<FMPosExt<positionClass>>& stack) override;

    /**
     * @brief Find the node path in the compressed De Bruijn graph corresponding
     * to a certain match
     *
     * @param occ The match of interest, containing its intervals in the SA and
     * reverse SA
     * @return FMOccSFI<positionClass> - An occurrence in the FM-index including
     * the node path in the compressed De Bruijn graph
     */
    FMOccSFI<positionClass>
    findNodePathForMatch(const FMOcc<positionClass>& occ);

    /**
     * @brief Find the node path in the compressed De Bruijn graph corresponding
     * to a certain match in a forward way
     *
     * @param pos The position corresponding to the match of interest
     * @param shift The shift corresponding to the match of interest
     * @param path The node path corresponding to the match of interest (to be
     * filled in)
     */
    void findNodePathForMatchForward(const positionClass& pos, const int& shift,
                                     std::vector<uint32_t>& path) const;

    /**
     * @brief Find to which strain a certain position in the text belongs
     *
     * @param input Position in the original text
     * @return int - Strain ID
     */
    int findStrain(length_t input);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR FILTERING STRAIN-FREE OCCURRENCES
    // ----------------------------------------------------------------------------

    /**
     * @brief Compare the edit distance of two strain-free occurrences in the
     * compressed de Bruijn graph.
     *
     * @param occ1 First occurrence to compare
     * @param occ2 Second occurrence to compare
     * @return true if the first occurrence has a smaller edit distance
     * @return false otherwise
     */
    static bool compareEditDistance(const FMOccSFR& occ1, const FMOccSFR& occ2);

    /**
     * @brief Alternative compare function used for sorting purposes. This
     * function first distinghuishes based on node path size.
     *
     * @param occ1 First occurrence to compare
     * @param occ2 Second occurrence to compare
     * @return true if the first occurrence is smaller than the second
     * @return false otherwise
     */
    static bool compareComplete(const FMOccSFR& occ1, const FMOccSFR& occ2);

    /**
     * @brief This function filters all occurrences that have the same node path
     * between themselves. When the function ends, the resulting set of
     * occurrences can contain maximum one occurrence corresponding to a certain
     * node path. Additionally, this function transforms the occurrences from
     * FMOcc<FMPosSFR> objects to FMOccSFR objects. This is necessary for the
     * further filtering process.
     *
     * @param occ The input vector of non-filtered occurrences
     * @param nonRedundantOcc The output vector for which the occurrences with
     * the same node path are filtered. This is to be filled in.
     * @param maxED The maximum edit distance
     * @param minLen The minimum length of the node paths of all occurrences
     * corresponding to this read. This is to be filled in.
     * @param maxLen The maximum length of the node paths of all occurrences
     * corresponding to this read. This is to be filled in.
     */
    void filterWithinSameNodePath(std::vector<FMOcc<FMPosSFR>>& occ,
                                  std::vector<FMOccSFR>& nonRedundantOcc,
                                  const int& maxED, size_t& minLen,
                                  size_t& maxLen) const;

    /**
     * @brief This function takes two occurrences and checks if the node path of
     * the shorter occurrence is a prefix of the node path of the longer
     * occurrence. It also reports the length of the longest common prefix of
     * the two node paths.
     *
     * @param prefixLength The length of the longest common prefix of the two
     * node paths. This is to be filled in.
     * @param shortestOcc The occurrence with the shortest node path, which
     * could be the prefix
     * @param longestOcc The occurrence with the longest node path
     * @return true the shortest node path is a complete prefix of the longest
     * one
     * @return false otherwise
     */
    bool checkNodePathsForPrefix(size_t& prefixLength, FMOccSFR* shortestOcc,
                                 FMOccSFR* longestOcc) const;

    /**
     * @brief This function clears the requested entries of vector
     * previousMatches by setting them to the nullpointer.
     *
     * @param previousMatches The vector of which the entries should be cleared
     * @param begin The begin index of the area to clear
     * @param end The end index of the area to clear (exclusive)
     */
    void clearPreviousMatches(std::vector<FMOccSFR*>& previousMatches,
                              size_t begin, size_t end) const;

    /**
     * @brief This function is called when a new minimum is found such that the
     * replacement related paramters for all occurrences in the current prefix
     * branch can be updated.
     *
     * @param previousMatches The vector containing the current prefix branch
     * @param previousMinimum The (new) minimum
     * @param previousMinimumChanged Boolean that indicates whether
     * previousMinimum has changed since the last call to this function.
     * @param index Index for the previousMatches vector indicating the
     * occurrence that is currently being considered.
     */
    void updateReplacements(std::vector<FMOccSFR*>& previousMatches,
                            FMOccSFR* previousMinimum,
                            bool previousMinimumChanged, size_t index) const;

    /**
     * @brief Update the value of previousMinimum to either shortestOcc or
     * longestOcc, depending on their edit distance.
     *
     * @param previousMinimumPointer Double pointer to the previous minimum. The
     * data of the first pointer or in other words, the second pointer will be
     * edited such that it points to the new minimum.
     * @param shortestOcc The candidate with the shortest node path
     * @param longestOcc The candidate with the longest node path
     */
    void updatePreviousMinimum(FMOccSFR** previousMinimumPointer,
                               FMOccSFR* shortestOcc,
                               FMOccSFR* longestOcc) const;

    /**
     * @brief This function is called when a shorter occurrence is found in the
     * prefix branch to compare the occurrence of interest with. It looks for a
     * prefix in the current prefix branch of the occurrence of interest. If a
     * prefix is found, the necessary values are updated.
     *
     * @param decrease The difference of the node path length of the occurrence
     * of interest and the occurrence in the prefix branch to consider first
     * @param m The current node of interest
     * @param previousMinimumPointer Double pointer to the previous minimum. The
     * data of the first pointer or in other words, the second pointer will be
     * edited such that it points to the new minimum if necessary.
     * @param previousMatches The vector containing the current prefix branch
     * @param previousMinimums The vector containing the old previous minimums
     * for all possible node path lengths
     * @param minLen The minimum length of the node paths of all occurrences
     * corresponding to this read
     * @param maxLen The maximum length of the node paths of all occurrences
     * corresponding to this read
     * @param len The length of the node path corresponding to the current
     * occurrence of interest
     * @param previousLen The length of the node path corresponding to the
     * previous occurrence of interest
     * @param maxED the maximum allowed edit distance
     * @return true if previousMinimum has changed
     * @return false otherwise.
     */
    bool handleIfPrefix(size_t& decrease, FMOccSFR* m,
                        FMOccSFR** previousMinimumPointer,
                        std::vector<FMOccSFR*>& previousMatches,
                        std::vector<FMOccSFR*>& previousMinimums,
                        const size_t& minLen, const size_t& maxLen,
                        const size_t& len, size_t& previousLen,
                        const int& maxED) const;

    /**
     * @brief This function considers one occurrence and checks if it has a
     * prefix that was previously considered. If so, it is analyzed whether it
     * replaces or is replaced by this other occurrence.
     *
     * @param occ The occurrence of interest
     * @param previousMinimumPointer Double pointer to the previous minimum. The
     * data of the first pointer or in other words, the second pointer will be
     * edited such that it points to the new minimum if necessary.
     * @param previousMatches The vector containing the current prefix branch
     * @param previousMinimums The vector containing the old previous minimums
     * for all possible node path lengths
     * @param minLen The minimum length of the node paths of all occurrences
     * corresponding to this read
     * @param maxLen The maximum length of the node paths of all occurrences
     * corresponding to this read
     * @param previousLen The length of the node path corresponding to the
     * previous occurrence of interest
     * @param maxED the maximum allowed edit distance
     */
    void handleReplacementConditionsForOccurrence(
        FMOccSFR& occ, FMOccSFR** previousMinimumPointer,
        std::vector<FMOccSFR*>& previousMatches,
        std::vector<FMOccSFR*>& previousMinimums, const size_t& minLen,
        const size_t& maxLen, size_t& previousLen, const int& maxED) const;

    /**
     * @brief This function iterates over all occurrences and analyzes them by
     * checking if the node paths are prefixes of eachother. If so, certain
     * attributes are set. The result of this function depends on the direction
     * of the node path (regular or reverse).
     *
     * @param nonRedundantOcc The set of not-filtered occurrences
     * @param minLen The minimum length of the node paths of all occurrences
     * corresponding to this read
     * @param maxLen The maximum length of the node paths of all occurrences
     * corresponding to this read
     * @param maxED the maximum allowed edit distance
     */
    void filterLinearInOneDirection(std::vector<FMOccSFR>& nonRedundantOcc,
                                    const size_t& minLen, const size_t& maxLen,
                                    const int& maxED) const;

    /**
     * @brief This function filters all occurrences in a linear manner. It is
     * more efficient than complete filtering, but does not eliminate all
     * redundant occurrences and rarely does not find the optimal replacement
     * for an occurrence. This function assumes that the occurrences were
     * already filtered such that for every possible node path, only the most
     * optimal occurrence was retained. This way, this function can compare
     * occurrences that have different node paths with each other.
     *
     * @param nonRedundantOcc The set of not-filtered occurrences
     * @param nonRedundantOccFinal The final set of filtered, non-redundant
     * occurrences. This is to be filled in.
     * @param maxED The maximum allowed edit distance
     * @param minLen The minimum length of the node paths of all occurrences
     * corresponding to this read
     * @param maxLen The maximum length of the node paths of all occurrences
     * corresponding to this read
     */
    void
    filterDifferentNodePathsLinear(std::vector<FMOccSFR>& nonRedundantOcc,
                                   std::vector<FMOccSFR>& nonRedundantOccFinal,
                                   const int& maxED, const size_t& minLen,
                                   const size_t& maxLen);

    /**
     * @brief This function filters all occurrences in a complete manner. It is
     * not very efficient, but eliminates all redundant occurrences and finds
     * the optimal replacement for every occurrence. This function assumes that
     * the occurrences were already filtered such that for every possible node
     * path, only the most optimal occurrence was retained. This way, this
     * function can compare occurrences that have different node paths with each
     * other.
     *
     * @param nonRedundantOcc The input list of occurrences that need to be
     * filtered. For every node path, only occurrence should be present.
     * @param nonRedundantOccFinal The output list of filtered occurrences. This
     * is to be filled in.
     * @param maxED The maximum allowed edit distance
     */
    void filterDifferentNodePathsComplete(
        std::vector<FMOccSFR>& nonRedundantOcc,
        std::vector<FMOccSFR>& nonRedundantOccFinal, const int& maxED) const;

    // ----------------------------------------------------------------------------
    // HELPER ROUTINES FOR VISUALIZATION
    // ----------------------------------------------------------------------------

    /**
     * @brief Initialize the outputfiles for the subgraph visualization process
     *
     * @param filename The base filename for all output files
     * @param multipleSubgraphs Indicates whether there are multiple subgraphs
     * corresponding to this pattern or not
     * @param edgefile The output stream corresponding to the edge file
     */
    void initializeFilesForVisualization(std::string filename,
                                         bool multipleSubgraphs,
                                         std::ofstream& edgefile);

    /**
     * @brief Make sure a node is visited in the visualization process
     *
     * @param id The node that needs to be visited
     * @param depth The current neigborhood depth
     * @param visited_nodes Array of nodes that have been visited
     * @param node_queue Queue of nodes that still need to be visited
     */
    void visitNode(uint32_t id, uint32_t depth,
                   std::vector<uint32_t>& visited_nodes,
                   std::queue<std::pair<int, int>>& node_queue);

    /**
     * @brief Create a visualization node and insert it in the vector of
     * visualization nodes
     *
     * @param visualizedNodes Vector containing all initialized visualization
     * nodes
     * @param id The node identifier of the node that needs to be visualized
     * @param node The node that needs to be visualized
     * @param path The original node path
     */
    void
    fillInVisualizationNode(std::vector<VisualizationNode*>& visualizedNodes,
                            uint32_t& id, Node& node,
                            std::vector<uint32_t>& path);

    /**
     * @brief Intermadiary function in the visualization process. It takes a
     * node and executes everything that is necessary to visualize this node
     * and its surroudings.
     *
     * @param path node path of the original match
     * @param subgraph_id ID of the subgraph (a match can have multiple
     * subgraphs)
     * @param visited_nodes array of nodes that have been visited
     * @param node_queue queue of nodes that still need to be visited
     * @param edgefile the output stream corresponding to the edge file
     * @param edgecounter counter for the key column
     * @param visualizedNodes Vector containing all initialized visualization
     * nodes
     */
    void visualizeSubgraphIntermediary(
        std::vector<uint32_t>& path, std::string subgraph_id,
        std::vector<uint32_t>& visited_nodes,
        std::queue<std::pair<int, int>>& node_queue, std::ofstream& edgefile,
        size_t& edgecounter, std::vector<VisualizationNode*>& visualizedNodes);

  public:
    // ----------------------------------------------------------------------------
    // ROUTINES FOR THE BUILDING AND LOADING PROCESS
    // ----------------------------------------------------------------------------

    /**
     * @brief Build the implicit representation of the compressed De Bruijn
     * graph based on the bidirectional FM-index
     *
     * @param baseFile base for the filenames
     * @param k parameter for the De Bruijn graph
     * @param sa_sparse sparseness factor for the SA
     * @param checkpoint_sparseness sparseness factors for the checkpoints that
     * aid in finding node identifiers.
     * @param progress prints progress if true
     */
    static void buildFMIndexDBG(const std::string& baseFile, const int k,
                                const std::vector<int>& sa_sparse,
                                const std::vector<int>& checkpoint_sparseness,
                                const bool progress) {
        FMIndexDBG(baseFile, k, sa_sparse, checkpoint_sparseness, progress);
    }

    /**
     * @brief Construct a new FMIndexDBG object from files in memory
     *
     * @param baseFile base for the filenames
     * @param sa_sparse sparseness factor for the SA
     * @param cp_sparse sparseness factor for the checkpoint k-mers
     * @param strainFree bool that indicates whether strain-free matching
     * will be used or not
     * @param option select algorithm option
     */
    FMIndexDBG(const std::string& baseFile, int sa_sparse, int cp_sparse,
               bool strainFree = false, bool filteringOptionComplete = false,
               const SelectOption& option = SelectOption::SIMPLE)
        : FMIndex<positionClass>(baseFile, sa_sparse) {
        std::cout << "Constructing the compressed de Bruijn graph..."
                  << std::endl;

        this->strainFree = strainFree;
        this->filteringOptionComplete = filteringOptionComplete;

        // Find the number of strains
        numberOfStrains =
            (std::count(this->bwt.begin(), this->bwt.end(), '%') + 1);
        this->numberOfSeparationCharacters = 2;

        // Find all start positions of the different strains and sort them along
        // with the position of the sentinel character
        sorted_startpositions.emplace_back(this->findSA(0));
        for (length_t i = 1; i < numberOfStrains; i++) {
            sorted_startpositions.emplace_back(this->findSA(i) + 1);
        }
        std::sort(sorted_startpositions.begin(), sorted_startpositions.end());

        std::string suffix =
            cp_sparse == INT32_MAX ? "none" : std::to_string(cp_sparse);

        // read the bit vectors for the implicit representation of the
        // compressed De Bruijn graph from files
        {
            std::ifstream ifs(baseFile + ".B.right." + suffix);
            if (!ifs) {
                throw std::runtime_error("Cannot open file " + baseFile +
                                         ".B.right." + suffix);
            }
            B_right.read(ifs);
            ifs.close();
        }
        {
            std::ifstream ifs(baseFile + ".B.left");
            if (!ifs) {
                throw std::runtime_error("Cannot open file " + baseFile +
                                         ".B.left");
            }
            B_left.read(ifs);
            ifs.close();
        }
        {
            std::ifstream ifs(baseFile + ".B.right.full." + suffix);
            if (!ifs) {
                throw std::runtime_error("Cannot open file " + baseFile +
                                         ".B.right.full." + suffix);
            }
            B_right_full.read(ifs);
            ifs.close();
        }
        // // Debugging
        // for (length_t i = 0; i < textLength; i++)
        // {
        //         std::cout << i << " B_left " << B_left[i]
        //         << std::endl;
        // }

        // read implicit representation of the compressed De Bruijn graph from
        // files
        {
            G = std::vector<Node>();
            Node n = Node();
            std::ifstream ifs(baseFile + ".DBG");
            if (!ifs) {
                throw std::runtime_error("Cannot open file " + baseFile +
                                         ".DBG");
            }
            ifs.read((char*)&k, sizeof(k));
            while (!n.read(ifs).eof()) {
                G.emplace_back(n);
            }
            ifs.close();
        }
        std::cout << "k = " << k << std::endl;
        k_DBG = k;
        numberOfGraphNodes = G.size();

        std::cout << "The pan-genome graph contains " << numberOfGraphNodes
                  << " nodes." << std::endl;

        // read the node identifier mappings from files
        {
            std::ifstream ifs(baseFile + ".right.map." + suffix,
                              std::ios::binary);
            if (!ifs) {
                throw std::runtime_error("Cannot open file " + baseFile +
                                         ".right.map." + suffix);
            }
            // ifs.seekg(0, std::ios::end);
            // mapping_right.resize(ifs.tellg() / sizeof(int)); //TODO can this
            // still be used?
            MappingPair mp = MappingPair();
            while (!mp.read(ifs).eof()) {
                mapping_right.emplace_back(mp);
            }
            ifs.close();
        }
        {
            std::ifstream ifs(baseFile + ".left.map", std::ios::binary);
            if (!ifs) {
                throw std::runtime_error("Cannot open file " + baseFile +
                                         ".left.map");
            }
            ifs.seekg(0, std::ios::end);
            mapping_left.resize(ifs.tellg() / sizeof(int));
            ifs.seekg(0, std::ios::beg);
            ifs.read((char*)&mapping_left[0],
                     mapping_left.size() * sizeof(int));
            ifs.close();
        }

        FMPos::setIndex(this);
    }

    // ----------------------------------------------------------------------------
    // GENERAL HELPER ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * @brief return whether we are matching strain-free or not
     *
     * @return true if we are matching strain-free
     * @return false otherwise
     */
    bool isStrainFree() const {
        return this->strainFree;
    }

    /**
     * @brief Return whether we are using complete filtering or not
     *
     * @return std::string - the filtering option
     */
    std::string getFilteringOption() const {
        if (filteringOptionComplete) {
            return "complete";
        } else {
            return "linear";
        }
    }

    /**
     * @brief Return the number of graph nodes that have been visited
     *
     * @return length_t - number of graph nodes visited
     */
    length_t getDBGNodes() const;

    /**
     * @brief Get the number of special cases in the linear filtering for
     * strain-free matching
     *
     * @return length_t - the number of special cases in the linear filtering
     * for strain-free matching
     */
    length_t getFilterSpecialCases() const;

    /**
     * @brief Resets the counters
     *
     */
    void resetCounters();

    /**
     * @brief Sets the search direction of the fm-index
     *
     * @param d the direction to search in, either FORWARD or BACKWARD
     */
    virtual void setDirection(Direction d) override {
        FMIndex<positionClass>::setDirection(d);
        FMPosSFR::setDirection(d);
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR MAPPING
    // ----------------------------------------------------------------------------

    /**
     * @brief If i is an index for the suffix array, find the node in the
     * compressed de Bruijn graph that corresponds to T[SA[i]..SA[i]+k[.
     * T[SA[j-1]..SA[j-1]+k[ can be at any location in the node.
     *
     * @param i Index in the suffix array indicating a certain k-mer of a
     * node
     * @param id the ID of the corresponding node in the compressed de Bruijn
     * graph (to be filled in)
     * @param l the number of characters that are before the k-mer in the node
     * (to be filled in)
     */
    void findID(length_t i, uint32_t& id, uint32_t& l) const;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR EXACT PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Exactly match the pattern to the compressed De Bruijn graph in a
     * strain-fixed way
     *
     * @param pattern the pattern to be matched
     * @return std::vector<TextOccurrenceSFI> - the result of the matching
     * procedure
     */
    std::vector<TextOccurrenceSFI> ExactMatchSFI(const std::string& pattern);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Approximately match the pattern to the compressed De Bruijn graph
     * using the naive backtracking strategy in a strain-fixed way
     *
     * @param pattern the pattern to be matched
     * @param maxED the maximum edit distance allowed
     * @return std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>> -
     * the results of the matching procedure
     */
    std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>
    approxMatchesNaiveSFI(const std::string& pattern, length_t maxED);

    /**
     * @brief This function maps matches in the SA along with their node path to
     * matches in the text. It takes the ranges of the matches and together with
     * the depth this is matched to a range in the text (this new range has a
     * width of depth). The edit distance is also mapped to this new range in
     * the text. This function also filters out the redundant matches using the
     * maximal allowed edit distance. This function is only used for
     * strain-fixed matching.
     *
     * @param occurrences the vector with all approximate matches and their
     * ranges in the SA and revSA
     * @param maxED the maximal allowed edit distance
     * @return std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>> -
     * a map of matches in the text containing a range and the edit distance as
     * values and the node path as key
     */
    std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>
    mapOccurrencesInSAToOccurrencesInTextSFI(
        std::vector<FMOcc<positionClass>>& occurrences, const int& maxED);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR FILTERING STRAIN-FREE OCCURRENCES
    // ----------------------------------------------------------------------------

    /**
     * @brief This function filters strain-free matches in the compressed de
     * Bruijn graph such that redundant matches are removed using the maximal
     * allowed edit distance.
     *
     * @param occurrences the vector with all strain-free approximate matches
     * @param maxED the maximal allowed edit distance
     * @return std::vector<FMOcc<FMPosSFR>> - the filtered set of occurrences
     */
    std::vector<FMOccSFR>
    filterStrainFreeMatches(std::vector<FMOcc<FMPosSFR>>& occurrences,
                            const int& maxED);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR VISUALIZATION
    // ----------------------------------------------------------------------------

    /**
     * @brief Display the subgraph corresponding to a certain path, given depth
     * > 0
     *
     * @param path the path that needs to visualized
     * @param depth the depth of the subgraph
     * @param filename the base for the filenames of the output files
     * @param multipleSubgraphs bool indicating whether there are multiple
     * subgraphs
     * @param subgraph_id id of the subgraph
     */
    void visualizeSubgraph(std::vector<uint32_t>& path, uint32_t depth,
                           std::string filename, bool multipleSubgraphs = false,
                           std::string subgraph_id = "");

    /**
     * @brief Visualize the subgraphs corresponding to the (approximate) matches
     * of a certain pattern. This function is used during strain-fixed
     * matching.
     *
     * @param paths map containing the node paths as keys and the corresponding
     * text occurrences as values
     * @param depth the depth of the subgraph
     * @param filename the base for the filenames of the output files
     */
    void visualizeSubgraphs(
        std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>& paths,
        uint32_t depth, std::string filename);

    /**
     * @brief Visualize the subgraphs corresponding to the (approximate) matches
     * of a certain pattern. This function is used during strain-free
     * matching.
     *
     * @param paths vector containing the strain-free matches
     * @param depth the depth of the subgraph
     * @param filename the base for the filenames of the output files
     */
    void visualizeSubgraphs(std::vector<FMOccSFR>& paths, uint32_t depth,
                            std::string filename);
};
