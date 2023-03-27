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

#include "fmindexDBG.h"
#include "divsufsort64.h"
#include "longestCommonPrefix.h"
#include "radix.h"
#include <set>
#include <sstream>

// ============================================================================
// CLASS FMIndexDBG
// ============================================================================

template <class positionClass>
thread_local length_t FMIndexDBG<positionClass>::nodeDBGCounter;
template <class positionClass>
thread_local length_t FMIndexDBG<positionClass>::filterSpecialCaseCounter;
template <class positionClass>
thread_local length_t FMIndex<positionClass>::nodeCounter;
template <class positionClass>
thread_local length_t FMIndex<positionClass>::matrixElementCounter;
template <class positionClass>
thread_local length_t FMIndex<positionClass>::positionsInPostProcessingCounter;
template <class positionClass>
thread_local length_t FMIndex<positionClass>::redundantNodePathsCounter;

// ----------------------------------------------------------------------------
// ROUTINES FOR THE BUILDING PROCESS
// ----------------------------------------------------------------------------

template <class positionClass>
void FMIndexDBG<positionClass>::createFMIndex(const std::string& baseFN,
                                              const std::vector<int>& sparse_sa,
                                              string& buf) {
    // read the text file from disk
    std::cout << "Reading " << baseFN << ".txt..." << std::endl;

    bool fastaInput = true;
    // If the buffer is empty, the index is not built from fasta files, but from
    // a preprocessed text file
    if (buf.empty()) {
        fastaInput = false;
        // Store the full text in a temporary buffer
        readTextOriginal(baseFN, buf);
    }

    // Find the length of the original text
    textLength = (buf[buf.size() - 1] == '\n') ? buf.size() - 1 : buf.size();

    bool newLine = false;

    // Remove return character if necessary
    if (textLength != buf.size()) {
        newLine = true;
        std::cout
            << "WARNING: the input text contained a tailing return, which was "
               "removed."
            << std::endl;
        buf = buf.substr(0, textLength);
    }

    // count the frequency of each characters in T
    std::vector<length_t> charCounts(256, 0);
    for (char c : buf)
        charCounts[(unsigned char)c]++;

    // count the number of unique characters in T
    int nUniqueChar = 0;
    for (length_t count : charCounts)
        if (count > 0)
            nUniqueChar++;

    std::cout << "\tText has length " << textLength << "\n";
    std::cout << "\tText has " << nUniqueChar << " unique characters\n";

    if (nUniqueChar > ALPHABET) {
        std::cerr << "FATAL ERROR: the number of unique characters in the "
                  << "text exceeds the alphabet size. Please recompile"
                  << "Nexus using a higher value for ALPHABET " << std::endl;
        exit(EXIT_FAILURE);
    }

    if (nUniqueChar < ALPHABET) {
        // Check if % is the character that is missing
        if (charCounts[37] == 0) {
            // The input should be a pan-genome with multiple strains for the
            // code to work.
            throw runtime_error(
                "Error: only one strain is present in the reference. As Nexus "
                "is tailored towards working with pan-genomes, matching "
                "patterns to one strain is currently not supported. If you "
                "really need this functionality though, then an empty strain "
                "should be added to resolve this issue.");
        }
        std::cout << "WARNING: the number of unique characters in the "
                  << "text is less than the ALPHABET size specified when "
                  << "Nexus was compiled. Performance may be affected\n";
    }

    // Create the alphabet
    sigma = Alphabet<ALPHABET>(charCounts);

    text = EncodedText<ALPHABET>(sigma, buf);
    text.write(baseFN + ".compressed.txt");

    std::cout << "Wrote file " << baseFN << ".compressed.txt\n";

    // write the character counts table
    {
        std::ofstream ofs(baseFN + ".cct", std::ios::binary);
        ofs.write((char*)charCounts.data(),
                  charCounts.size() * sizeof(length_t));
        ofs.close();
    }

    std::cout << "Wrote file " << baseFN << ".cct\n";

    // Create cumulative character counts
    length_t cumCount = 0;
    for (size_t i = 0; i < charCounts.size(); i++) {
        if (charCounts[i] == 0)
            continue;
        counts.push_back(cumCount);
        cumCount += charCounts[i];
    }

    // build the SA
    std::cout << "Generating the suffix array using divsufsort64...\n";
    clock_t startTime = clock();
    int64_t* SA = (int64_t*)malloc((size_t)textLength * sizeof(int64_t));
    divsufsort64((uchar*)&(buf)[0], SA, textLength);
    clock_t endTime = clock();
    float milliseconds = clock_diff_to_msec(endTime - startTime);
    printf("divsufsort64 took [%.2fs]\n", milliseconds / 1000.0);

    // perform a sanity check on the suffix array
    std::cout << "\tPerforming sanity checks..." << std::endl;
    // sanityCheck(text, SA);
    // std::cout << "\tSanity checks OK" << std::endl;

    // If the index is built from fasta files, then we don't clear the buffer
    // because otherwise all of the preprocessing should be redone. So building
    // the index from a preprocessed text is somewhat more memory-efficient
    // right now.
    if (!fastaInput) {
        // Clear the buffer to save space
        buf.clear();
        buf.resize(0);
        buf.shrink_to_fit();
    }

    // build the BWT
    std::cout << "Generating BWT..." << std::endl;
    bwt.resize(textLength);
    for (size_t i = 0; i < textLength; i++)
        if (SA[i] > 0)
            bwt.set(i, text[SA[i] - 1]);
        else
            bwt.set(i, text[textLength - 1]);

    // Write bwt
    bwt.write(baseFN + ".bwt");
    std::cout << "Wrote file " << baseFN << ".bwt\n";

    // create succinct BWT bitvector table
    fwdRepr = BWTRepr<ALPHABET>(sigma, bwt);
    fwdRepr.write(baseFN + ".brt");
    std::cout << "Wrote file: " << baseFN << ".brt" << std::endl;

    // Count the number of strains
    numberOfStrains = fwdRepr.occ(sigma.c2i('%'), textLength) + 1;
    bwt.clear();
    fwdRepr.clear();

    // create sparse suffix arrays
    for (int saSF : sparse_sa) {
        sparseSA.clear();
        sparseSA = SparseSuffixArray(SA, saSF, textLength);
        sparseSA.write(baseFN);
        std::cout << "Wrote sparse suffix array with factor " << saSF
                  << std::endl;
    }
    delete[] SA;

    sparseSA.clear();

    // Reverse the original text
    std::cout << "Reversing the original text...\n";
    startTime = clock();

    // If we are not working with fasta files, then the buffer was cleared to
    // save space. So the text should be read from memory again.
    if (!fastaInput) {
        readTextOriginal(baseFN, buf);

        if (newLine) {
            buf.pop_back();
        }
    }

    reverse(buf.begin(), buf.end());
    endTime = clock();
    milliseconds = clock_diff_to_msec(endTime - startTime);
    printf("Reversing took [%.2fs]\n", milliseconds / 1000.0);

    // build the reverse SA
    std::cout << "Generating the reverse suffix array using divsufsort64...\n";
    startTime = clock();
    int64_t* revSA = (int64_t*)malloc((size_t)textLength * sizeof(int64_t));
    divsufsort64((uchar*)&(buf)[0], revSA, textLength);
    endTime = clock();
    milliseconds = clock_diff_to_msec(endTime - startTime);
    printf("divsufsort64 took [%.2fs]\n", milliseconds / 1000.0);

    buf.clear();
    buf.resize(0);
    buf.shrink_to_fit();

    // perform a sanity check on the suffix array
    std::cout << "\tPerforming sanity checks..." << std::endl;
    // sanityCheck(text, revSA);
    // std::cout << "\tSanity checks OK" << std::endl;

    // build the reverse BWT
    std::cout << "Generating reverse BWT..." << std::endl;
    revbwt.resize(textLength);
    for (size_t i = 0; i < textLength; i++)
        if (revSA[i] > 0)
            revbwt.set(i, text[textLength - revSA[i]]);
        else
            revbwt.set(i, text[0]);

    // Encode reverse bwt
    revbwt.write(baseFN + ".rev.bwt");
    delete[] revSA;

    std::cout << "Wrote file " << baseFN << ".rev.bwt\n";

    // create succinct reverse BWT bitvector table
    revRepr = BWTRepr<ALPHABET>(sigma, revbwt);
    revRepr.write(baseFN + ".rev.brt");
    std::cout << "Wrote file: " << baseFN << ".rev.brt" << std::endl;

    revbwt.clear();
    revRepr.clear();
    text.clear();
}

template <class positionClass>
void FMIndexDBG<positionClass>::computeLCPKasai(Bitvec2& LCP, bool progress,
                                                uint& k) {
    std::cout << "Computing the LCP..." << std::endl;
    // Initialize the LCP array with one more element than the length of the
    // original text
    LCP = (textLength + 1);
    // Create a rank array, which will contain the inverse of the suffix array
    std::vector<length_t> rank(textLength);
    // Fill in the rank array
    for (length_t i = 0; i < textLength; i++) {
        rank[findSA(i)] = i;
        if (progress) {
            if (i % (textLength / 100) == 0) {
                std::cout << "Progress part 1/5: " << i / (textLength / 100)
                          << "%"
                          << "\r";
                std::cout.flush();
            }
        }
    }
    if (progress) {
        std::cout << "Progress part 1/5: " << 100 << "%"
                  << "\n";
    }
    uint h = 0;
    // Fill in the LCP array by iterating over all suffixes in sorted order
    // (this is the order in the suffix array)
    for (length_t i = 0; i < textLength; i++) {
        // Verify that the suffix starting at index i is not the
        // lexicographically smallest suffix. In practice, the lexicographically
        // smallest suffix is "$", so this one is skipped.
        if (rank[i] > 0) {
            // The suffix starting at index j is successor of the suffix
            // starting at index i in the suffix array.
            length_t j = findSA(rank[i] - 1);
            // Increment h until it represents the length of the longest common
            // prefix between the two suffixes
            while (text[i + h] == text[j + h]) {
                h++;
            }
            // Update the compacted LCP array.
            if (h < k) {
                LCP[rank[i]] = 0;
            } else if (h == k) {
                LCP[rank[i]] = 1;
            } else {
                LCP[rank[i]] = 2;
            }
            // Decrement h if necessary. This is to avoid redundant computations
            // in the next iteration.
            if (h > 0) {
                h--;
            }
        }
        // Print progress
        if (progress) {
            if (i % (textLength / 100) == 0) {
                std::cout << "Progress part 2/5: " << i / (textLength / 100)
                          << "%"
                          << "\r";
                std::cout.flush();
            }
        }
    }
    // Print progress
    if (progress) {
        std::cout << "Progress part 2/5: " << 100 << "%"
                  << "\n";
    }
    // for (length_t i = 0; i < textLength + 1; i++)
    // {
    //         std::cout << i << " " << LCP[i] << std::endl;
    // }
}

template <class positionClass>
void FMIndexDBG<positionClass>::computeBitVectors(std::queue<int>& Q,
                                                  Bitvec& Br_right,
                                                  Bitvec& Bl_right,
                                                  bool progress, uint& k) {
    //  Build the longest common prefix array
    Bitvec2 LCP;
    computeLCPPrezza(this, LCP, progress, k);

    // Build the bit vectors
    std::cout << "Computing Br_right and Bl_right..." << std::endl;

    // Initialize all values to false
    for (length_t i = 0; i < textLength; i++) {
        Br_right[i] = false;
        Bl_right[i] = false;
    }

    // Create a set of all positions in the SA that are within a distance of k
    // of a separation character. This is necessary for the correct encoding of
    // end nodes.
    std::set<length_t> end_values;
    for (length_t i = 0; i < numberOfStrains; i++) {
        length_t index = i;
        end_values.insert(index);
        for (uint j = 1; j < k; j++) {
            index = findLF(index, false);
            end_values.insert(index);
        }
    }
    // lb stores the last index for which the LCP value is lower then k
    length_t lb = 0;
    // kIndex stores the last index for which the LCP value is equal to k
    length_t kIndex = 0;
    // lastDiff stores the last index at which the characters BWT[lastDiff-1]
    // and BWT[lastDiff] differ
    length_t lastDiff = 0;
    // if open is true, we are in an area for which the LCP values are larger
    // than k
    bool open = false;
    // a counter to create node identifiers
    int counter = 0;
    // create a copy of the cumulative counts array corresponding to the
    // bidirectional FM-index
    std::vector<length_t> counts_copy = counts;
    // Iterate over all entries in the LCP array
    for (length_t i = 1; i < textLength + 1; i++) {
        // the counts array is adjusted every iteration so that it would contain
        // the following information: counts[BWT[j]] =LF(j)+1, where j is the
        // index of the last occurrence of character BWT[j] in BWT before the
        // current index i of the algorithm
        counts_copy[bwt[i - 1]]++;

        if (LCP[i] >= 1) {
            // We enter a range with LCP entries bigger than or equal to k
            open = true;
            if (LCP[i] == 1) {
                // update the last index for which the LCP value is equal to k
                kIndex = i;
            }
        } else {
            if (open) {
                // We are leaving a k-mer range that must be analyzed

                // Check if the k-mer is right-maximal
                if (kIndex > lb) {
                    // LCP[j] = k for at least one j with j bigger than the last
                    // index for which the LCP value is lower then k

                    // We are now in a right-maximal k-mer interval

                    // First we check that the current position in the SA does
                    // not correspond to an end node
                    if (end_values.find(lb) == end_values.end()) {
                        // We are not situated in an end node

                        // We set the boundaries of the closed SA interval
                        // corresponding to the right-maximal k-mer to true in
                        // Br_right
                        Br_right[lb] = true;
                        Br_right[i - 1] = true;
                        // We get the range in the reverse SA that corresponds
                        // to this k-mer
                        Range reverseRange = getReverseRangeKMer(lb, false, k);
                        length_t lb_rev = reverseRange.getBegin();
                        // Add the new node to the graph
                        G.emplace_back(k, i - lb, lb, lb, lb_rev, lb_rev);
                        // Push the node identifier to the queue
                        Q.push(counter);
                        // Increment the node ID creator
                        counter++;
                    }
                }
                // check whether the k-mer corresponding to the current interval
                // is a left-maximal repeat by checking if the last index
                // lastDiff at which the characters BWT[lastDiff-1] and
                // BWT[lastDiff] differ, is higher than the left boundary (lb)
                // of the SA interval
                if (lastDiff > lb) {
                    // Iterate over all characters in the BWT corresponding to
                    // this left-maximal k-mer interval
                    for (length_t j = lb; j < i; j++) {
                        // Get the corresponding character

                        char c = sigma.i2c(bwt[j]);
                        // Check if the preceding character is a separation
                        // character, which means we are situated in a start
                        // node
                        if (c != '$' && c != '%') {
                            // This is not a start node
                            int cIdx = sigma.c2i(c);
                            // the counts array can directly be used to find the
                            // right boundary of the closed SA interval of the
                            // preceding k-mer, which needs to be set to true in
                            // Bl_right
                            length_t index = counts_copy[cIdx] - 1;
                            // Also check that we are not situated in an end
                            // node
                            if (end_values.find(index) == end_values.end()) {
                                // We are not in an end node, so the bit can be
                                // set to true
                                Bl_right[index] = true;
                            }
                        }
                    }
                }
                // We have left the k-mer interval
                open = false;
            }
            // update the last index for which the LCP value is lower then k
            lb = i;
        }
        if (bwt[i] != bwt[i - 1]) {
            // the last index at which the characters BWT[lastDiff-1] and
            // BWT[lastDiff] differ
            lastDiff = i;
        }
        // Print progress
        if (progress) {
            if (i % (textLength / 100) == 0) {
                std::cout << "Progress part 1/3: " << i / (textLength / 100)
                          << "%"
                          << "\r";
                std::cout.flush();
            }
        }
    }
    // Print progress
    if (progress) {
        std::cout << "Progress part 1/3: " << 100 << "%"
                  << "\n";
    }
    open = false;
    // the one-bits in Bl_right that also correspond to a right-maximal k-mer
    // must be reset to zero
    for (length_t i = 0; i < textLength; i++) {
        if (open) {
            // We are in a right-maximal k-mer interval
            Bl_right[i] = false;
            if (Br_right[i]) {
                // We exit a right-maximal k-mer interval
                open = false;
            }
        } else if (Br_right[i]) {
            // We enter a right-maximal k-mer interval
            Bl_right[i] = false;
            open = true;
        }
        // Print progress
        if (progress) {
            if (i % (textLength / 100) == 0) {
                std::cout << "Progress part 2/3: " << i / (textLength / 100)
                          << "%"
                          << "\r";
                std::cout.flush();
            }
        }
    }

    // Index the bit vectors such that they support rank operations
    Br_right.index();
    Bl_right.index();

    // Print progress
    if (progress) {
        std::cout << "Progress part 2/3: " << 100 << "%"
                  << "\n";
    }

    // // Debugging
    // for (length_t i = 0; i < textLength; i++) {
    //     std::cout << i << " Br_right " << Br_right[i] << " Bl_right "
    //               << Bl_right[i] << std::endl;
    // }
}

template <class positionClass>
std::vector<std::pair<char, Range>>
FMIndexDBG<positionClass>::getIntervals(length_t i, length_t j) {

    std::vector<std::pair<char, Range>> result;
    result.reserve(sigma.size());
    // Iterate over all characters
    for (length_t k = 0; k < sigma.size(); k++) {
        SARangePair child = SARangePair(Range(i, j), Range(0, 0));
        // Find the range over the SA after appending the extra character to the
        // front of the current suffix
        if (findRangesWithExtraCharBackward(
                k, SARangePair(Range(i, j), Range(0, 0)), child)) {
            result.emplace_back(sigma.i2c(k), child.getRangeSA());
        }
    }
    return result;
}

template <class positionClass>
Range FMIndexDBG<positionClass>::getReverseRangeKMer(length_t indexInSA,
                                                     bool isEndNode, uint& k) {
    // Find the index of the substring of interest in the original text
    length_t indexInText = findSA(indexInSA);
    // Get the substring of interest: the first k-mer of the suffix at the
    // indexInSA'th position in the SA. Match this k-mer to the reference
    SARangePair p = matchStringBidirectionally(Substring(
        (text.decodeSubstring(sigma, indexInText, (indexInText + k)))));
    // If the node is not an end node, all computations are done
    if (!isEndNode) {
        // Return the reverse SA range
        return p.getRangeSARev();
    }
    // If the node is an end node, multiple end nodes can exist with the same
    // corresponding substring. Hence, we need to distinguish between these as
    // well. Following variable keeps track of the number of extra characters we
    // use to distinguish the nodes.
    length_t extra = textLength / 100 + k;
    // Keep adding sets of characters until the match is unique or the start of
    // the reference was reached
    while (p.width() > 1 && extra <= indexInText) {
        // Exactly match the substring with extra characters
        p = matchStringBidirectionally(Substring(
            text.decodeSubstring(sigma, indexInText - extra, indexInText + k)));
        extra += textLength / 100 + k;
    }
    // If too much characters were added in the last set, match again until the
    // start of the reference
    if (extra > indexInText) {
        p = matchStringBidirectionally(
            Substring(text.decodeSubstring(sigma, 0, indexInText + k)));
    }
    // Return the first entry of the found reverse SA range. We don't return the
    // complete reverse SA range because it could contain more than one entry if
    // for example the first two strains are identical.
    return Range(p.getRangeSARev().getBegin(),
                 p.getRangeSARev().getBegin() + 1);
}

template <class positionClass>
void FMIndexDBG<positionClass>::buildCompressedGraph(
    const std::vector<int>& checkpoint_sparseness, bool progress,
    std::vector<RankInterface>& B_rights,
    std::vector<std::vector<MappingPair>>& mapping_rights, uint& k) {
    std::cout << "Start building the implicit graph:" << std::endl;

    // Counter for progress printing
    int counter = 0;
    // Initialize the Br_right and Bl_right bit vectors to the length of the
    // text + 1. This way, rank(textLength) can be called as well. These are
    // temporary bit vectors for the build process.
    Bitvec Br_right = (textLength + 1);
    Bitvec Bl_right = (textLength + 1);

    const int cp_len = checkpoint_sparseness.size();

    // Set the start values of B_right and B_left to 0. These are
    // the bit vectors that will be stored.
    for (length_t i = 0; i < textLength; i++) {
        B_right[i] = false;
        B_left[i] = false;
        for (int j = 0; j < cp_len; j++) {
            B_rights[j][i] = false;
        }
    }
    // queue Q stores which nodes still need to be handled by the construction
    // algorithm
    std::queue<int> Q;
    // Store a temporary mapping of node IDs with index of the leftmost k-mer
    // in the reverse SA
    std::vector<pair<length_t, int>> temporary_mapping_left;
    // Store temporary mappings of node IDs and offsets with index of the
    // corresponding k-mer in the SA
    std::vector<std::vector<std::pair<length_t, MappingPair>>>
        temporary_mapping_rights(checkpoint_sparseness.size());

    // Compute temporary bit vectors Br_right and Bl_right
    computeBitVectors(Q, Br_right, Bl_right, progress, k);

    // Number of right-maximal nodes
    int rightMax = Br_right.rank(textLength) / 2;
    // Number of nodes that are not right-maximal and precede a left-maximal
    // node
    int leftMax = Bl_right.rank(textLength);

    G.reserve(rightMax + leftMax + numberOfStrains);

    std::cout << "Constructing the compressed graph..." << std::endl;

    // Add dummy nodes that will represent the nodes that are not right-maximal
    // and precede a left-maximal node later
    for (int s = 0; s < leftMax; s++) {
        G.emplace_back();
    }

    // Add the end nodes (they can still be extended later)
    for (length_t s = 0; s < numberOfStrains; s++) {
        int id = rightMax + leftMax + s;
        Range reverseRange = getReverseRangeKMer(s, true, k);
        G.emplace_back(1, 1, s, s, reverseRange.getBegin(), 0);
        Q.push(id);
        Bl_right[s] = 0;
    }

    // Find the number of graph nodes
    numberOfGraphNodes = G.size();
    // Forward matching
    setDirection(FORWARD);

    // Keep iterating over the queue until all nodes are handled
    while (!Q.empty()) {
        // Get the next node from the queue
        int id = Q.front();
        Node node = G[id];
        Q.pop();

        // extendable keeps track of whether the node can still be extended to
        // the left (true) or the end of the node was reached (false)
        bool extendable = true;
        // lb represents the left bound of the SA ranges of the total string
        // corresponding to the current node
        length_t lb = node.left_kmer_forward;

        // Keep extending the node to the left until the beginning of the node
        // is reached
        while (extendable) {
            extendable = false;

            // Get all possible intervals that result from extending the current
            // interval with one character to the left
            std::vector<std::pair<char, Range>> list =
                getIntervals(lb, lb + node.multiplicity);
            // Iterate over all resulting intervals along with the corresponding
            // character
            for (std::pair<char, Range> p : list) {
                // Extract the character that was used to extend the match
                char c = p.first;
                // Get the bounds of the new SA range
                length_t i = p.second.getBegin();
                length_t j = p.second.getEnd();
                // Check whether the k-length prefix of the suffix at index i in
                // the SA is a right-maximal k-mer. If so, this node is complete
                // and we can process the next interval or node
                int ones = Br_right.rank(i + 1);
                if (ones % 2 == 0 && Br_right[i] == 0) {
                    // the k-length prefix of the suffix at index i in the SA is
                    // a not right-maximal k-mer

                    // Check that the next character is not a separation
                    // character, since these should not occur in the middle of
                    // a node
                    if (c != '%' && c != '$') {
                        if (list.size() == 1) {
                            // Only one character could be used for extension,
                            // so the node can be extended
                            extendable = true;
                            node.len = node.len + 1;
                            node.left_kmer_forward = i;
                            lb = i;

                            // Add checkpoint if necessary
                            if (node.len > k) {
                                for (int j = 0; j < cp_len; j++) {
                                    if ((node.len - k) %
                                            checkpoint_sparseness[j] ==
                                        0) {
                                        // Set the correct bit in B_right
                                        B_rights[j][node.left_kmer_forward +
                                                    node.multiplicity - 1] =
                                            true;

                                        // Update the temporary right mapping
                                        temporary_mapping_rights[j]
                                            .emplace_back(
                                                node.left_kmer_forward +
                                                    node.multiplicity - 1,
                                                MappingPair(id, node.len - k));
                                    }
                                }
                            }

                        } else {
                            // Multiple characters can be used for extension, so
                            // the node is left-maximal. This means we have
                            // reached the beginning of the node.

                            // The predecessor of this node needs to be
                            // initialized
                            int newID = rightMax + Bl_right.rank(i);
                            Range reverseRange = getReverseRangeKMer(
                                i, newID >= rightMax + leftMax, k);
                            G[newID] =
                                Node(k, j - i, i, i, reverseRange.getBegin(),
                                     reverseRange.getBegin());
                            Q.push(newID);
                        }
                    }
                }
            }
        }

        // If the node was extended or the node is an end node, index of the
        // leftmost k-mer in the reverse suffix array must be updated.
        if (node.len > k || id >= rightMax + leftMax) {
            Range reverseRange = getReverseRangeKMer(
                node.left_kmer_forward, id >= rightMax + leftMax, k);
            node.left_kmer_reverse = reverseRange.getBegin();
        }

        // Set the correct bit in B_left corresponding to this node
        length_t indexForBit = node.left_kmer_reverse + node.multiplicity - 1;
        B_left[indexForBit] = true;
        // Update the temporary left mapping
        temporary_mapping_left.emplace_back(indexForBit, id);

        // Set the correct bit in B_right corresponding to this node
        indexForBit = node.right_kmer_forward + node.multiplicity - 1;
        for (int j = 0; j < cp_len; j++) {
            B_rights[j][indexForBit] = true;
            // Update the temporary right mapping
            temporary_mapping_rights[j].emplace_back(indexForBit,
                                                     MappingPair(id, 0));
        }

        if (node.len > k) {
            // Set the correct bit in B_right corresponding to this node
            indexForBit = node.left_kmer_forward + node.multiplicity - 1;
            for (int j = 0; j < cp_len; j++) {

                if ((node.len - k) % checkpoint_sparseness[j] > 0) {
                    B_rights[j][indexForBit] = true;
                    // Update the temporary right mapping
                    temporary_mapping_rights[j].emplace_back(
                        indexForBit, MappingPair(id, node.len - k));
                }
            }
        }

        // Update the node attributes in the graph
        G[id] = node;

        // Print progress
        if (progress && numberOfGraphNodes > 100) {
            if (counter % (numberOfGraphNodes / 100) == 0) {
                std::cout << "Progress part 3/3: "
                          << counter / (numberOfGraphNodes / 100) << "%"
                          << "\r";
                std::cout.flush();
            }
        }
        counter++;
    }
    // Print progress
    if (progress) {
        std::cout << "Progress part 3/3: " << 100 << "%"
                  << "\n";
    }

    // Index the three vectors such that they support rank operations. B_right
    // also supports select operations.
    B_left.indexInterface();
    for (int j = 0; j < cp_len; j++) {
        B_rights[j].indexInterface();
    }

    // Sanity check
    if (size_t(rightMax + leftMax + numberOfStrains) !=
        B_left.rank(textLength)) {
        std::cout << "Warning: node count is incorrect." << std::endl;
    }

    // // Debugging
    // for (length_t i = 0; i < textLength; i++) {
    //     std::cout << i << " B_right " << B_rights[0][i] << " B_left "
    //               << B_left[i] << std::endl;
    // }

    // Sanity check
    if (numberOfGraphNodes != temporary_mapping_left.size()) {
        std::cout << "Warning: node count is incorrect." << std::endl;
    }

    // Build the mapping that maps identifiers deduced from B_left to the true
    // identifiers
    mapping_left.resize(numberOfGraphNodes);
    for (pair<length_t, int> p : temporary_mapping_left) {
        int id_left = B_left.rank(p.first);
        mapping_left[id_left] = p.second;
    }

    for (int j = 0; j < cp_len; j++) {
        // Build the mapping that maps identifiers deduced from B_right to the
        // true identifiers
        mapping_rights[j].resize(temporary_mapping_rights[j].size());
        for (pair<length_t, MappingPair> p : temporary_mapping_rights[j]) {
            int id_right = B_rights[j].rank(p.first);
            mapping_rights[j][id_right] = p.second;
        }
    }
}

// ----------------------------------------------------------------------------
// HELPER ROUTINES FOR MAPPING
// ----------------------------------------------------------------------------

template <class positionClass>
length_t FMIndexDBG<positionClass>::getDBGNodes() const {
    return nodeDBGCounter;
}

template <class positionClass>
length_t FMIndexDBG<positionClass>::getFilterSpecialCases() const {
    return filterSpecialCaseCounter;
}

template <class positionClass> void FMIndexDBG<positionClass>::resetCounters() {
    nodeCounter = 0;
    matrixElementCounter = 0;
    positionsInPostProcessingCounter = 0;
    redundantNodePathsCounter = 0;
    nodeDBGCounter = 0;
    filterSpecialCaseCounter = 0;
    elapsedNodePaths = std::chrono::duration<double>::zero();
    elapsedSAtoText = std::chrono::duration<double>::zero();
}

template <class positionClass>
int FMIndexDBG<positionClass>::findIDLast(length_t i) const {
    // Find the number of ones up to i in B_right
    int id = B_right.rank(i);
    nodeDBGCounter++;
    // Retrieve the correct node ID using the node mapping
    return mapping_right[id].id;
}

template <class positionClass>
int FMIndexDBG<positionClass>::findIDFirst(length_t i) const {
    // Find the number of ones up to i in B_left
    int id = B_left.rank(i);
    nodeDBGCounter++;
    // Retrieve the correct node ID using the node mapping
    return mapping_left[id];
}

template <class positionClass>
void FMIndexDBG<positionClass>::findID(length_t j, uint32_t& id,
                                       uint32_t& l) const {
    nodeDBGCounter++;
    j--;
    l = 0;
    // The offset of the checkpoint in the node
    uint32_t l_offset, id_right;
    // Move to the left until an identifier is found
    while (B_right[j] == 0) {
        // Shift the k-length window one character to the left
        j = findLF(j, false);
        l++;
    }
    // Find the number of ones up to i in B_right
    id_right = B_right.rank(j);
    // Retrieve the correct node ID using the node mapping
    id = mapping_right[id_right].id;
    // Retrieve the offset of the checkpoint in the node
    l_offset = mapping_right[id_right].distanceFromRightEnd;
    l = G[id].len - k_DBG - l_offset + l;
    return;
}

template <class positionClass>
int FMIndexDBG<positionClass>::jumpToSuccessorThroughEdge(
    uint32_t id, uint32_t& offset_reverse) const {
    // Move one step to the right
    length_t left_kmer_reverse =
        findLF(G[id].right_kmer_reverse + offset_reverse, true);
    // Find the ID of the successor
    id = findIDFirst(left_kmer_reverse);
    // Find the new reverse offset
    offset_reverse = left_kmer_reverse - G[id].left_kmer_reverse;
    return id;
}

template <class positionClass>
int FMIndexDBG<positionClass>::jumpToPredecessorThroughEdge(
    uint32_t id, uint32_t& offset) const {
    // Move one step to the right
    length_t right_kmer_forward =
        findLF(G[id].left_kmer_forward + offset, false);
    // Find the ID of the successor
    id = findIDLast(right_kmer_forward);
    // Find the new offset
    offset = right_kmer_forward - G[id].right_kmer_forward;
    return id;
}

template <class positionClass>
bool FMIndexDBG<positionClass>::jumpToSuccessorWithChar(
    uint32_t id, uint32_t& id_successor, uint32_t posInAlphabet,
    uint32_t reverse_offset) const {
    Node node = G[id];
    // Create a range pair that represents the substring corresponding to the
    // current node in the graph
    SARangePair p(Range(node.left_kmer_forward,
                        node.left_kmer_forward + node.multiplicity),
                  Range(node.right_kmer_reverse,
                        node.right_kmer_reverse + node.multiplicity));
    // Append the next character to the substring and find the corresponding
    // range pair
    if (!findRangesWithExtraCharForward(posInAlphabet, p, p)) {
        // No successor was found
        return false;
    }
    // offset is only used in end nodes ending with character '%'. It checks
    // that the offset is not larger than the number of such end nodes that can
    // be accessed.

    if (reverse_offset >= p.width()) {
        // No extra successor was found
        return false;
    }
    // Find the identifier of the successor we are looking for.
    id_successor = findIDFirst(p.getRangeSARev().getBegin() + reverse_offset);
    return true;
}

template <class positionClass>
bool FMIndexDBG<positionClass>::jumpToPredecessorWithChar(
    uint32_t id, uint32_t& id_predecessor, uint32_t posInAlphabet,
    uint32_t offset) const {
    Node node = G[id];
    // Create a range pair that represents the substring corresponding to the
    // current node in the graph
    SARangePair p(Range(node.left_kmer_forward,
                        node.left_kmer_forward + node.multiplicity),
                  Range(node.right_kmer_reverse,
                        node.right_kmer_reverse + node.multiplicity));
    // Prepend the next character to the substring and find the corresponding
    // range pair
    if (!findRangesWithExtraCharBackward(posInAlphabet, p, p)) {
        // No predecessor was found
        return false;
    }
    // Find the identifier of the predecessor we are looking for.
    id_predecessor = findIDLast(p.getRangeSA().getBegin());
    return true;
}

template <class positionClass>
std::vector<TextOccurrenceSFI>
FMIndexDBG<positionClass>::convertToMatchesInTextSFI(
    const FMOccSFI<positionClass>& saMatch) {

    return convertToMatchesInTextSFI(saMatch.getRanges(), saMatch.getNodePath(),
                                     saMatch.getDistanceFromLeftEnd(),
                                     saMatch.getDepth(), saMatch.getDistance(),
                                     saMatch.getShift());
}

template <class positionClass>
std::vector<TextOccurrenceSFI>
FMIndexDBG<positionClass>::convertToMatchesInTextSFI(
    const SARangePair& ranges, const std::vector<uint32_t>& nodepath,
    const uint32_t& distanceFromLeftEnd, const int& patternLength,
    const int& distance, const length_t& shift) {
    std::vector<TextOccurrenceSFI> textMatches;
    textMatches.reserve(ranges.width());

    for (length_t i = ranges.getRangeSA().getBegin();
         i < ranges.getRangeSA().getEnd(); i++) {
        // find the startPosition in the text by looking at the SA
        length_t startPos = findSA(i) + shift;

        // cap startPos at textLength
        startPos = startPos % textLength;

        length_t endPos = startPos + patternLength;

        textMatches.emplace_back(Range(startPos, endPos), distance, nodepath,
                                 distanceFromLeftEnd);
    }
    return textMatches;
}

template <class positionClass>
const bool
FMIndexDBG<positionClass>::separationIsNext(positionClass pos) const {
    return pos.separationIsNext();
}

template <class positionClass>
void FMIndexDBG<positionClass>::extendFMPos(
    const SARangePair& parentRanges,
    std::vector<FMPosExt<positionClass>>& stack, int row, int trueDepth) const {

    // iterate over the entire alphabet, excluding the separation characters
    for (length_t i = 2; i < sigma.size(); i++) {

        extendFMPosIntermediary(parentRanges, stack, row, i, trueDepth);
    }
}

template <class positionClass>
void FMIndexDBG<positionClass>::extendFMPos(
    const positionClass& pos,
    std::vector<FMPosExt<positionClass>>& stack) const {

    pos.extendFMPos(stack);
}

template <class positionClass>
FMOccSFI<positionClass> FMIndexDBG<positionClass>::findNodePathForMatch(
    const FMOcc<positionClass>& occ) {
    auto start = chrono::high_resolution_clock::now();
    // Refer to the correct path finding function
    std::vector<uint32_t> path = {};
    uint32_t distanceFromLeftEnd = 0;
    findNodePathForMatchForward(occ.getPosition(), occ.getShift(), path,
                                distanceFromLeftEnd);
    redundantNodePathsCounter++;
    auto finish = chrono::high_resolution_clock::now();
    elapsedNodePaths += finish - start;
    return FMOccSFI<positionClass>(occ, path, distanceFromLeftEnd);
}

template <class positionClass>
void FMIndexDBG<positionClass>::findNodeUnderK(
    const positionClass& pos, const int& shift,
    std::vector<uint32_t>&
        path) { // TODO: greedy version, could use some more tweaking
    length_t originalDepth = pos.getDepth();

    setDirection(BACKWARD);

    std::vector<FMPosExt<positionClass>> stack;
    std::vector<FMPosExt<positionClass>> candidatesAfterSearch;
    std::vector<FMPosExt<positionClass>> candidatesNearStartOfGenome;
    stack.push_back(FMPosExt<positionClass>((char)0, pos));
    length_t originalWidth = pos.getRanges().width();
    while (!stack.empty()) {
        FMPosExt<positionClass> p = stack.back();
        stack.pop_back();
        if (p.getTrueDepth() == k_DBG) {
            candidatesAfterSearch.push_back(p);
            originalWidth -= p.getRanges().width();
        } else {
            extendFMPos(p, stack);
        }
        length_t newWidth = 0;
        for (auto tempPos : stack) {
            newWidth += tempPos.getRanges().width();
        }
        if (newWidth != originalWidth) {
            candidatesNearStartOfGenome.push_back(p);
        }
    }

    setDirection(FORWARD);

    for (FMPosExt<positionClass> p : candidatesNearStartOfGenome) {
        stack.push_back(FMPosExt<positionClass>((char)0, p));
    }

    while (!stack.empty()) {
        FMPosExt<positionClass> p = stack.back();
        stack.pop_back();
        if (p.getTrueDepth() == k_DBG) {
            candidatesAfterSearch.push_back(p);
            originalWidth -= p.getRanges().width();
        } else {
            extendFMPos(p, stack);
        }
        length_t newWidth = 0;
        for (auto tempPos : stack) {
            newWidth += tempPos.getRanges().width();
        }
        if (newWidth != originalWidth) {
            candidatesNearStartOfGenome.push_back(p);
        }
    }

    std::vector<pair<uint32_t, uint32_t>> stack2;

    for (FMPosExt<positionClass> p : candidatesAfterSearch) {
        // Get the left boundary of the SA interval of the match
        length_t rb = p.getRanges().getRangeSA().getEnd();
        // Find the corresponding node of the path along with the
        // distance of the k-mer of the match to the beginning of the node
        uint32_t id, l;
        findID(rb, id, l);
        if (!std::count(path.begin(), path.end(), id)) {
            path.push_back(id);
            auto len = G[id].len;
            auto var = len - l - k_DBG + originalDepth;
            if (var < k_DBG) {
                stack2.push_back(make_pair(id, var - originalDepth));
            }
        }
    }

    while (!stack2.empty()) {
        auto pair = stack2.back();
        stack2.pop_back();
        for (size_t i = 2; i < sigma.size(); i++) {
            // Check if there exists a successor that is the result of
            // extension with the current character
            uint32_t id_successor;
            if (jumpToSuccessorWithChar(pair.first, id_successor, i)) {
                if (!std::count(path.begin(), path.end(), id_successor)) {
                    path.push_back(id_successor);
                    auto var = pair.second + G[id_successor].len - k_DBG + 1 +
                               originalDepth;
                    if (var < k_DBG) {
                        stack2.push_back(
                            make_pair(id_successor, var - originalDepth));
                    }
                }
            }
        }
    }
}

template <class positionClass>
void FMIndexDBG<positionClass>::findNodePathForMatchForward(
    const positionClass& occ, const int& shift, std::vector<uint32_t>& path,
    uint32_t& distanceFromLeftEnd) {
    // A match shorter than k does not have a unique node path
    if (occ.getTrueDepth() < k_DBG) {
        // Find all corresponding nodes
        findNodeUnderK(occ, shift, path);
        return;
    }
    // Get the left boundary of the SA interval of the match
    length_t lb = occ.getRanges().getRangeSA().getBegin();
    // Get the start position of an occurrence in the reference text
    length_t positionInText = findSA(lb) + shift;
    // Get the string corresponding to the first k-mer of the occurrence
    string firstKmer =
        text.decodeSubstring(sigma, positionInText, positionInText + k_DBG);
    // Get the SA range corresponding to this first k-mer
    Range kmerRange = matchString(firstKmer);
    // Find the first node of the path along with the distance of the first
    // k-mer of the match to the beginning of the first node
    uint32_t id, pos;
    findID(kmerRange.getEnd(), id, pos);
    // Set distanceFromLeftEnd
    distanceFromLeftEnd = pos;
    // Find the distance of the the beginning of first k-mer of the match to the
    // end of the first node. This reflects how much of the occurrence has been
    // handled.
    pos = G[id].len - pos;
    // Store the start node in the node path
    path.emplace_back(id);
    while (pos < occ.getDepth()) {
        // Copy the old id
        uint32_t oldID = id;
        // Find the ID of the successor using the next character of the
        // occurrence
        jumpToSuccessorWithChar(oldID, id, text[positionInText + pos]);
        // Find the new position in the match based on the length of the new
        // node
        pos += G[id].len - (k_DBG - 1);
        // Add the new node to the path
        path.emplace_back(id);
    }
}

template <class positionClass>
int FMIndexDBG<positionClass>::findStrain(length_t input) const {
    int b = 0;
    int e = sorted_startpositions.size();
    if (input >= sorted_startpositions[e - 1]) {
        throw std::runtime_error(
            "Binary search for strain: this input is not allowed");
    }
    // Take an educated first guess
    length_t avgLength = textLength / numberOfStrains;
    int i = input / avgLength;
    // Use binary search to find the correct strain
    while (b < e) {
        if (sorted_startpositions[i] <= input &&
            sorted_startpositions[i + 1] > input) {
            return i;
        } else if (sorted_startpositions[i] > input) {
            e = i;
        } else {
            b = i;
        }
        i = floor((b + e) / 2);
    }
    throw std::runtime_error("Could not finish binary search for strain");
}

// ----------------------------------------------------------------------------
// ROUTINES FOR MAPPING
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

template <class positionClass>
std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>
FMIndexDBG<positionClass>::approxMatchesNaiveSFI(const std::string& pattern,
                                                 length_t maxED) {

    std::vector<FMOcc<positionClass>> occurrences =
        approxMatchesNaiveIntermediate(pattern, maxED);

    return mapOccurrencesInSAToOccurrencesInTextSFI(occurrences, maxED);
}

template <class positionClass>
std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>
FMIndexDBG<positionClass>::mapOccurrencesInSAToOccurrencesInTextSFI(
    std::vector<FMOcc<positionClass>>& occ, const int& maxED) {
    auto start = chrono::high_resolution_clock::now();

    sort(occ.begin(), occ.end());
    occ.erase(unique(occ.begin(), occ.end()), occ.end());
    std::vector<TextOccurrenceSFI> occurrencesInText;
    occurrencesInText.reserve(1000 * maxED);
    // map the startposition to the best match (lowest edit distance)
    map<length_t, TextOccurrenceSFI> posToBestMatch;

    if (occ.size() == 0) {
        auto finish = chrono::high_resolution_clock::now();
        elapsedSAtoText += finish - start;
        return {};
    }
    if (occ.size() == 1) {
        // all occ are distinct
        positionsInPostProcessingCounter = occ[0].getRanges().width();
        const auto& m = convertToMatchesInTextSFI(findNodePathForMatch(occ[0]));
        std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>> paths;
        for (const auto& mOcc : m) {
            paths[mOcc.getNodePath()].emplace_back(mOcc);
            auto& newOcc = paths[mOcc.getNodePath()].back();
            newOcc.generateOutput();
            newOcc.setStrain(findStrain(newOcc.getRange().getBegin()));
        }
        for (const auto& myPair : paths) {
            sort(paths[myPair.first].begin(), paths[myPair.first].end());
        }
        auto finish = chrono::high_resolution_clock::now();
        elapsedSAtoText += finish - start;
        return paths;
    }

    // more than 1 reported occurrence, could be redundant
    for (const auto& it : occ) {

        const Range& range = it.getRanges().getRangeSA();
        positionsInPostProcessingCounter += range.width();

        const auto& matchesInTextToCheck =
            convertToMatchesInTextSFI(findNodePathForMatch(it));
        occurrencesInText.insert(occurrencesInText.end(),
                                 matchesInTextToCheck.begin(),
                                 matchesInTextToCheck.end());
    }

    sort(occurrencesInText.begin(), occurrencesInText.end());

    std::vector<TextOccurrenceSFI> nonRedundantOcc;
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
            // check if this later occurrence is better than the previous one
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

    std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>> paths;
    for (TextOccurrenceSFI& occ : nonRedundantOcc) {
        occ.generateOutput();
        occ.setStrain(findStrain(occ.getRange().getBegin()));
        paths[occ.getNodePath()].emplace_back(occ);
    }
    auto finish = chrono::high_resolution_clock::now();
    elapsedSAtoText += finish - start;
    return paths;
}

template <class positionClass>
bool FMIndexDBG<positionClass>::compareEditDistance(const FMOccSFR& occ1,
                                                    const FMOccSFR& occ2) {
    return occ1.getDistance() < occ2.getDistance();
}

template <class positionClass>
bool FMIndexDBG<positionClass>::compareComplete(const FMOccSFR& occ1,
                                                const FMOccSFR& occ2) {
    const auto& pos1 = occ1.getPosition();
    const auto& path1 = pos1.getNodePath();
    const auto& pos2 = occ2.getPosition();
    const auto& path2 = pos2.getNodePath();
    if (path1.size() != path2.size()) {
        // The node paths do not have the same size, shorter is smaller
        return path1.size() < path2.size();
    } else if (path1 != path2) {
        // The node paths are not identical, they are compared elementwise to
        // see which is smaller
        return path1 < path2;
    } else if (pos1.getDistanceFromLeftEnd() != pos2.getDistanceFromLeftEnd()) {
        // The distance of the left end of the matches to the left end of the
        // first node of the nodepath is not identical, smaller distance is
        // smaller in general
        return pos1.getDistanceFromLeftEnd() < pos2.getDistanceFromLeftEnd();
    } else if (occ1.getDistance() != occ2.getDistance()) {
        // begin is equal, better ed is smarter
        return occ1.getDistance() < occ2.getDistance();
    } else if (pos1.getTrueDepth() != pos2.getTrueDepth()) {
        // begin is equal, shorter match is smaller
        return pos1.getTrueDepth() < pos2.getTrueDepth();
    } else {
        // prefer no shift
        return occ1.getShift() < occ2.getShift();
    }
}

template <class positionClass>
void FMIndexDBG<positionClass>::filterWithinSameNodePath(
    std::vector<FMOcc<FMPosSFR>>& occ, std::vector<FMOccSFR>& nonRedundantOcc,
    const int& maxED, size_t& minLen, size_t& maxLen) const {

    length_t maxDiff = 2 * maxED;
    length_t prevBegin = numeric_limits<length_t>::max();
    int prevED = maxED + 1;
    length_t prevDepth = numeric_limits<length_t>::max();

    const std::vector<uint32_t>* prevNodePath = nullptr;

    // initialize the variables that store the minimum and maximum node path
    // sizes for this occurrence.
    minLen = occ[0].getPosition().getNodePath().size();
    maxLen = occ[0].getPosition().getNodePath().size();

    // Iterate over all occurrences to remove redundant ones with the same node
    // path
    for (const auto& o : occ) {
        size_t len = o.getPosition().getNodePath().size();
        if (len > maxLen) {
            // update the maximum node path length
            maxLen = len;
        }
        if (len < minLen) {
            // update the minimum node path length
            minLen = len;
        }
        if (prevNodePath && (*prevNodePath != o.getPosition().getNodePath())) {
            // An occurrence with a new node path is always added

            prevBegin = o.getPosition().getDistanceFromLeftEnd();
            prevED = o.getDistance();
            prevDepth = o.getPosition().getTrueDepth();
            prevNodePath = &o.getPosition().getNodePath();

            nonRedundantOcc.emplace_back(o);
            continue;
        }
        auto diff = abs_diff<length_t>(o.getPosition().getDistanceFromLeftEnd(),
                                       prevBegin);
        if (diff == 0) {
            continue;
        }
        if (diff <= maxDiff) {
            // check if this later occurrence is better than the previous one
            if (o.getDistance() > prevED) {
                continue;
            }
            if (o.getDistance() == prevED &&
                o.getPosition().getTrueDepth() >= prevDepth) {
                continue;
            }

            // prev was worse so pop_back
            nonRedundantOcc.pop_back();
        }

        prevBegin = o.getPosition().getDistanceFromLeftEnd();
        prevED = o.getDistance();
        prevDepth = o.getPosition().getTrueDepth();
        prevNodePath = &o.getPosition().getNodePath();

        nonRedundantOcc.emplace_back(o);
    }
}

template <class positionClass>
bool FMIndexDBG<positionClass>::checkNodePathsForPrefix(
    size_t& prefixLength, FMOccSFR* shortestOcc, FMOccSFR* longestOcc) const {
    // extract the node path of the shortest occurrence
    const std::vector<uint32_t>& shorterNodepath =
        shortestOcc->getPosition().getNodePath();
    // extract the node path of the longest occurrence
    const std::vector<uint32_t>& longerNodepath =
        longestOcc->getPosition().getNodePath();
    // Check if the shortest node path is a prefix of the longer one
    for (prefixLength = 0; prefixLength < shorterNodepath.size();
         prefixLength++) {
        // Every element is checked
        if (shorterNodepath[prefixLength] != longerNodepath[prefixLength]) {
            return false;
        }
    }
    // All checks passed, the shortest path is a prefix
    return true;
}

template <class positionClass>
void FMIndexDBG<positionClass>::clearPreviousMatches(
    std::vector<FMOccSFR*>& previousMatches, size_t begin, size_t end) const {
    // Set the correct entries to null pointers
    std::fill(previousMatches.begin() + begin, previousMatches.begin() + end,
              nullptr);
}

template <class positionClass>
void FMIndexDBG<positionClass>::updateReplacements(
    std::vector<FMOccSFR*>& previousMatches, FMOccSFR* previousMinimum,
    bool previousMinimumChanged, size_t index) const {
    if (!previousMinimumChanged) {
        // In this case, previousMinimum did not change since the previous call
        // to updateReplacement. Hence, only the occurrence that was just
        // handled needs to be updated, not all occurrences in the prefixbranch.

        // Check if previousMinimum is a better replacement for
        // previousMatches[index] and if so, update the replacement attributes
        // for previousMatches[index]
        if ((FMOccSFR::setReplacement)(previousMatches[index],
                                       previousMinimum)) {
            // If previousMinimum is a better replacement for
            // previousMatches[index], indicate that previousMinimum replaces
            // previousMatches[index]
            (FMOccSFR::replaces)(previousMinimum, previousMatches[index]);
        }
        return;
    }

    // In this case, previousMinimum did change since the previous call to
    // updateReplacement. Hence, all occurrences in the current prefixbranch
    // must be updated.

    // Iterate over all occurrences that are currently being observed
    for (size_t i = 0; i < previousMatches.size(); i++) {
        // Check that both pointers are valid
        if (previousMatches[i] && previousMinimum) {
            // Check if previousMinimum is a better replacement for
            // previousMatches[i] and if so, update the replacement attributes
            // for previousMatches[i]
            if ((FMOccSFR::setReplacement)(previousMatches[i],
                                           previousMinimum)) {
                // If previousMinimum is a better replacement for
                // previousMatches[i], indicate that previousMinimum replaces
                // previousMatches[i]
                (FMOccSFR::replaces)(previousMinimum, previousMatches[i]);
            }
        }
    }
}

template <class positionClass>
void FMIndexDBG<positionClass>::updatePreviousMinimum(
    FMOccSFR** previousMinimumPointer, FMOccSFR* shortestOcc,
    FMOccSFR* longestOcc) const {

    if (longestOcc->getDistance() < shortestOcc->getDistance()) {
        // The longer occurrence is better
        *previousMinimumPointer = longestOcc;
    } else {
        // The longer occurrence has a larger or equal edit distance. Hence, the
        // shorter occurrence is better
        *previousMinimumPointer = shortestOcc;
    }
}

template <class positionClass>
bool FMIndexDBG<positionClass>::handleIfPrefix(
    size_t& decrease, FMOccSFR* m, FMOccSFR** previousMinimumPointer,
    std::vector<FMOccSFR*>& previousMatches,
    std::vector<FMOccSFR*>& previousMinimums, const size_t& minLen,
    const size_t& maxLen, const size_t& len, size_t& previousLen,
    const int& maxED) const {

    // Extract pointer to the previousMinimum FMOccSFR object
    FMOccSFR* previousMinimum = *previousMinimumPointer;

    // Following variable will store the length of the longest common prefix of
    // the two node paths of interest
    size_t prefixLength;
    // Check if the node path of previousMatches[len - minLen - decrease] is a
    // prefix of that of occurrence m
    if (checkNodePathsForPrefix(prefixLength,
                                previousMatches[len - minLen - decrease], m)) {
        // A complete prefix was found

        // Now check if one of these two occurrences can actually replace
        // each other. This is only the case if their starting point wrt the
        // first node is within a maximum distance.
        if (!(FMOccSFR::checkProximity)(
                previousMatches[len - minLen - decrease], m, maxED)) {
            // The occurrences are too far apart, and cannot replace each other

            // Count the amount of times this special case occurs
            filterSpecialCaseCounter++;
            // Since the new occurrence does not belong to the current prefix
            // branch, a new one must be started. Therefore, the previousMinimum
            // must be reported to represent the
            // previous branch.
            (FMOccSFR::report)(previousMinimum);
            // All entries before m are cleared as we have established that they
            // are not part of the current prefix branch.
            clearPreviousMatches(previousMatches, 0, len - minLen);
            // As there is no prefix, m is the new previous minimum.
            *previousMinimumPointer = m;
            // previousMinimum was changed.
            return true;
        }

        if (len <= previousLen) {
            // The previous occurrence that was considered had an equal or
            // larger node length wrt this occurrence. This means that we are
            // not continuing in the same line of prefixes, but we are creating
            // a new branch. For this reason, the previous minimum must be
            // reported to represent the previous branch.
            (FMOccSFR::report)(previousMinimum);
        }

        // We check which of the two occurrences currently being considered is
        // better and update the previous minimum accordingly.
        updatePreviousMinimum(previousMinimumPointer,
                              previousMinimums[len - decrease - minLen], m);
        // Return true if previousMinimum was changed.
        return *previousMinimumPointer != previousMinimum;
    }

    // In this case, the node path of previousMinimums[len - decrease - minLen]
    // is not a prefix of that of m.

    // Since the node path is not a complete prefix, a new branch of prefixes is
    // created. Hence the previousMinimum must be reported to represent the
    // previous branch.
    (FMOccSFR::report)(previousMinimum);

    // The variable below is necessary to call checkNodePathsForPrefix but
    // otherwise not of interest
    size_t i;
    // If the current read also has occurrences with a node path of length
    // prefixLength (as previously found), we check if the current prefixbranch
    // contains an entry that has this node path length. If so, we check if the
    // node path corresponding to this entry is a prefix of the node path of m.
    if (prefixLength >= minLen && previousMatches[prefixLength - minLen] &&
        checkNodePathsForPrefix(i, previousMatches[prefixLength - minLen], m)) {
        // We have found an shorter entry in the current prefixbranch that is a
        // complete prefix of occurrence m.

        // We clear all entries that have a node path length between
        // prefixLength and len (the length of the node path of m), since we
        // have established that these are not prefixes of m and hence do not
        // belong in the current prefix branch.
        clearPreviousMatches(previousMatches, prefixLength - minLen + 1,
                             len - minLen);
        // We check which of the two occurrences currently being considered is
        // better and update the previous minimum accordingly.
        updatePreviousMinimum(previousMinimumPointer,
                              previousMinimums[prefixLength - minLen], m);
        // Return true if previousMinimum was changed.
        return *previousMinimumPointer != previousMinimum;
    }

    // In this case, no entry in previousMatches was found that is a prefix of
    // m. For this reason, all entries before m are cleared as we have
    // established that they are not part of the current prefix branch.
    clearPreviousMatches(previousMatches, 0, len - minLen);
    // As there is no prefix, m is the new previous minimum.
    *previousMinimumPointer = m;
    // previousMinimum was changed.
    return true;
}

template <class positionClass>
void FMIndexDBG<positionClass>::handleReplacementConditionsForOccurrence(
    FMOccSFR& occ, FMOccSFR** previousMinimumPointer,
    std::vector<FMOccSFR*>& previousMatches,
    std::vector<FMOccSFR*>& previousMinimums, const size_t& minLen,
    const size_t& maxLen, size_t& previousLen, const int& maxED) const {

    // Extract pointer to the previousMinimum FMOccSFR object
    FMOccSFR* previousMinimum = *previousMinimumPointer;

    // Extract pointer to the current FMOccSFR object of interest
    FMOccSFR* m = &occ;
    // Store the length of the node path corresponding to the current FMOccSFR
    // object of interest
    size_t len = occ.getPosition().getNodePath().size();

    // Store the current FMOccSFR object of interest in the previousMatches to
    // add it to the prefix branch
    previousMatches[len - minLen] = m;
    // Clear all entries in previousMatches that have a longer node path than m,
    // because these are not part of the current prefix branch. Otherwise, they
    // would have been considered after m due to the sorting definition.
    clearPreviousMatches(previousMatches, len - minLen + 1,
                         maxLen - minLen + 1);

    // The variable below indicates which shorter entry in previousMatches we
    // are considering to be compared with m. Decrease indicates the difference
    // of the node path length between the two occurrences.
    size_t decrease = 1;
    // The boolean below indicates whether a shorter entry in previousMatches
    // was found to compare m with.
    bool otherMatchFound = false;
    // The boolean below indicates whether the previousMinimum pointer has
    // changed or not.
    bool previousMinimumChanged = true;
    // Iterate previousMatches from right to left (longer to shorter node paths)
    // starting at the entry before the entry corresponding to m. Stop iterating
    // if the beginning of previousMatches is reached.
    while (len - decrease >= minLen) {
        // Check if the entry in previousMatches of interest is a valid pointer
        if (!previousMatches[len - minLen - decrease]) {
            // No valid pointer, so the decrease value is increased such that
            // the previous entry can be checked in the next iteration
            decrease++;
            // Continue the while loop to keep searching for a shorter
            // occurrence to compare this occurrence with.
            continue;
        }

        // In this case, a valid entry in previousMatches is found. Therefore
        // otherMatchFound is set to true.
        otherMatchFound = true;
        // For the current entry in previousMatches, handleIfPrefix looks for a
        // prefix in the current prefix branch and handles things accordingly.
        // previousMinimumChanged is also updated.
        previousMinimumChanged = handleIfPrefix(
            decrease, m, previousMinimumPointer, previousMatches,
            previousMinimums, minLen, maxLen, len, previousLen, maxED);
        // Stop the while loop here, since a valid entry was already found.
        break;
    }

    if (!otherMatchFound) {
        // If no valid entry was found during the while loop, there was no
        // shorter occurrence to compare with. Hence, we are starting a new
        // prefixbranch. Therefore, the previous minimum must be reported such
        // that the previous branch is represented.
        (FMOccSFR::report)(previousMinimum);
        // Since there is no other occurrence for comparison, m is the new
        // previous minimum.
        *previousMinimumPointer = m;
    }

    // Store the new minimum in the correct place in previousMinimums (based on
    // its node path length), such that it can later be used again.
    previousMinimums[len - minLen] = *previousMinimumPointer;
    // Update the occurrences in the current prefixbranch as wel as the
    // attributes of the minimum of this branch such that replacement
    // information is stored.
    updateReplacements(previousMatches, *previousMinimumPointer,
                       previousMinimumChanged, len - minLen);
    // Update the value that stores the length of the previously considered node
    // path.
    previousLen = len;
}

template <class positionClass>
void FMIndexDBG<positionClass>::filterLinearInOneDirection(
    std::vector<FMOccSFR>& nonRedundantOcc, const size_t& minLen,
    const size_t& maxLen, const int& maxED) const {

    // Create the previousMatches vector, which will store the prefixbranch that
    // is considered at a given time in the algorithm. A prefix branch is a set
    // of occurrences, all of which are prefixes of each other.
    std::vector<FMOccSFR*> previousMatches(maxLen - minLen + 1, nullptr);
    // Create the previousMinimums vector, which will store the minimum
    // occurrences that correspond to the occurrences in previousMatches. A
    // minimum occurrence in previousMinimums is always has a shorter or equal
    // node path length wrt the corresponding entry in previousMatches.
    std::vector<FMOccSFR*> previousMinimums(maxLen - minLen + 1, nullptr);

    // Create a dummy occurrences to serve as the first previousMinimum. This
    // way, the previousMinimum pointer is always valid, so no checks on that
    // need to be present.
    FMOccSFR nullOcc = FMOccSFR();
    // Create the previousMinimum pointer
    FMOccSFR* previousMinimum = &nullOcc;

    // Initialize the variable that stores the length of the occurrence that was
    // considered previously.
    size_t previousLen = 0;

    // Iterate over all occurrences corresponding to the current read of
    // interest.
    for (FMOccSFR& occ : nonRedundantOcc) {
        // For the occurrence that is currently considered, analyze if it
        // replaces or is replaced by another occurrence that was previously
        // considered.
        handleReplacementConditionsForOccurrence(
            occ, &previousMinimum, previousMatches, previousMinimums, minLen,
            maxLen, previousLen, maxED);
    }

    // Report the final previousMinimum, such that the last prefix branch is
    // also represented.
    (FMOccSFR::report)(previousMinimum);
}

template <class positionClass>
void FMIndexDBG<positionClass>::filterDifferentNodePathsLinear(
    std::vector<FMOccSFR>& nonRedundantOcc,
    std::vector<FMOccSFR>& nonRedundantOccFinal, const int& maxED,
    const size_t& minLen, const size_t& maxLen) {

    // The number of iterations over all occurrences is currently hard coded.
    // This should be either 2 or 3. The final value is still to be chosen.
    const int nrOfIterations = 3;

    // Start the first filtering iteration on the regular node paths.
    filterLinearInOneDirection(nonRedundantOcc, minLen, maxLen, maxED);

    // Change the direction in which the node paths should be analyzed to
    // reverse
    FMOccSFR::setDirection(true);

    // Reverse all node paths and reset the correct attributes
    for (FMOccSFR& m : nonRedundantOcc) {
        (FMOccSFR::reverseNodePath)(m);
    }

    // Sort the occurrences again, since the node paths are reversed
    sort(nonRedundantOcc.begin(), nonRedundantOcc.end());

    // Start the second filtering iteration on the reverse node paths.
    filterLinearInOneDirection(nonRedundantOcc, minLen, maxLen, maxED);

    // Change the direction in which the node paths should be analyzed to
    // regular again
    FMOccSFR::setDirection(false);

    if (nrOfIterations == 3) {

        // Reverse all node paths and reset the correct attributes
        for (FMOccSFR& m : nonRedundantOcc) {
            (FMOccSFR::reverseNodePath)(m);
        }

        // Sort the occurrences again, since the node paths are reversed back to
        // regular
        sort(nonRedundantOcc.begin(), nonRedundantOcc.end());

        // Start the third filtering iteration on the regular node paths. This
        // time, we have extra information from the second filtering iteration.
        filterLinearInOneDirection(nonRedundantOcc, minLen, maxLen, maxED);
    }

    // Iterate over all occurrences and add only the ones that are reported as
    // non-redundant by the linear sorting algorithm to the final list of
    // occurrences.
    for (FMOccSFR& m : nonRedundantOcc) {
        if (m.isNonRedundant()) {
            nonRedundantOccFinal.emplace_back(m);
        }
    }

    if (nrOfIterations == 2) {

        // Reverse all reported node paths (and reset the correct attributes)
        for (FMOccSFR& m : nonRedundantOccFinal) {
            (FMOccSFR::reverseNodePath)(m);
        }
    }
}

template <class positionClass>
void FMIndexDBG<positionClass>::filterDifferentNodePathsComplete(
    std::vector<FMOccSFR>& nonRedundantOcc,
    std::vector<FMOccSFR>& nonRedundantOccFinal, const int& maxED) const {

    // Iterate over all allowed edit distances, starting at 0. The reason for
    // this is occurrences with a smaller edit distance are always the first
    // choice as replacement for for other occurrences.
    for (int ed = 0; ed <= maxED; ed++) {
        // Now iterate over all remaining pairs of occurrences to remove
        // redundant occurrences with a different node path using two for-loops

        // The for loop to select the first occurrence of the pair. This will be
        // the replacing occurrence if applicable.
        for (size_t i = 0; i < nonRedundantOcc.size(); i++) {
            // Only consider this replacing occurrence if its edit distance is
            // equal to the one that is considered by the outer loop.
            if (nonRedundantOcc[i].getDistance() == ed) {
                // Store the node path of the replacing occurrence.
                const std::vector<uint32_t>* nodepath_i =
                    &nonRedundantOcc[i].getPosition().getNodePath();
                // The for loop to select the second occurrence of the pair.
                // This will be the replaced occurrence if applicable (due to
                // how the occurrences were sorted earlier).
                for (size_t j = 0; j < nonRedundantOcc.size(); j++) {
                    // Store the node path of the replaced occurrence.
                    const std::vector<uint32_t>* nodepath_j =
                        &nonRedundantOcc[j].getPosition().getNodePath();
                    // Check that the replaced occurrence was not already
                    // replaced by another occurrence earlier. If so, we do not
                    // need to continue searching for a worse replacement (due
                    // to sorting). The node path lengths must also be
                    // different, otherwise it is not possible that one node
                    // path is part of the other. Next, check that it either has
                    // a larger edit distance than the replacing occurrence, or
                    // the same edit distance but a longer node path (due to
                    // sorting).
                    if (!nonRedundantOcc[j].isToBeRemoved() &&
                        nodepath_i->size() != nodepath_j->size() &&
                        ((nonRedundantOcc[j].getDistance() > ed) ||
                         (nonRedundantOcc[j].getDistance() >= ed && j > i))) {
                        // It is still possible that nonRedundantOcc[i] replaces
                        // nonRedundantOcc[j].

                        // Create two new node path pointers
                        const std::vector<uint32_t>*nodepath1, *nodepath2;
                        // Store the shortest node path in nodepath1, the
                        // longest node path in nodepath2
                        if (nodepath_i->size() < nodepath_j->size()) {
                            nodepath1 = nodepath_i;
                            nodepath2 = nodepath_j;
                        } else {
                            nodepath1 = nodepath_j;
                            nodepath2 = nodepath_i;
                        }

                        // Redundant occurrences have a different node path
                        // length at this point. Node path 1 has the shortest
                        // length.

                        // Find the last possible index in node path 2 where
                        // node path 1 could start
                        std::vector<uint32_t>::const_iterator lastIndex =
                            nodepath2->end() - nodepath1->size();

                        // Check if the first element of node path 1 occurs in
                        // node path 2 and if so, where
                        std::vector<uint32_t>::const_iterator iter =
                            std::find_if(nodepath2->begin(), lastIndex + 1,
                                         [nodepath1](const uint32_t id) {
                                             return id == (*nodepath1)[0];
                                         });

                        // Check if node path 1 could still fit in node path 2
                        // starting from the index that was found
                        while (iter < lastIndex + 1) {

                            // Find the index relative to the start of node path
                            // 2
                            unsigned int index =
                                std::distance(nodepath2->begin(), iter);

                            // The node at index iter in node path 2 is equal to
                            // the first node of node path 1
                            bool equal = true;

                            // Check if node path 1 occurs at this location
                            // completely
                            for (unsigned int l = 1; l < nodepath1->size();
                                 l++) {
                                if ((*nodepath1)[l] !=
                                    (*nodepath2)[index + l]) {
                                    // No occurrence of node path 1
                                    equal = false;
                                    break;
                                }
                            }

                            // One node path effectively occurs within the other
                            if (equal) {
                                // One node path effectively occurs within the
                                // other

                                // Now check if one of these two occurrences can
                                // actually replace each other. This is only the
                                // case if their starting point wrt the graph
                                // is within a maximum distance.

                                // Offset 1 stores the distance of the start of
                                // the shortest occurrence to the start of the
                                // longest path
                                // Offset 2 stores the distance of the start of
                                // the longest occurrence to the start of the
                                // corresponding longest path
                                uint32_t offset1, offset2;
                                if (nodepath_i->size() < nodepath_j->size()) {
                                    offset1 = nonRedundantOcc[i]
                                                  .getPosition()
                                                  .getDistanceFromLeftEnd();
                                    offset2 = nonRedundantOcc[j]
                                                  .getPosition()
                                                  .getDistanceFromLeftEnd();
                                } else {
                                    offset1 = nonRedundantOcc[j]
                                                  .getPosition()
                                                  .getDistanceFromLeftEnd();
                                    offset2 = nonRedundantOcc[i]
                                                  .getPosition()
                                                  .getDistanceFromLeftEnd();
                                }
                                // Also take the nodes of the longest path that
                                // occurs before the shortest path in
                                // consideration
                                for (size_t nodeIterator = 0;
                                     nodeIterator < index; nodeIterator++) {
                                    offset1 +=
                                        G[(*nodepath2)[nodeIterator]].len -
                                        k_DBG + 1;
                                }

                                length_t maxDiff = 2 * maxED;

                                auto diff =
                                    abs_diff<length_t>(offset1, offset2);

                                // Check if the distance between the occurrences
                                // is within the maximum allowed distance
                                if (diff <= maxDiff) {
                                    // nonRedundantOcc[j] can be replaced by
                                    // nonRedundantOcc[i]
                                    nonRedundantOcc[j].setToBeRemoved();
                                    // nonRedundantOcc[i] replaces
                                    // nonRedundantOcc[j]
                                    nonRedundantOcc[i].setToBeKept();
                                }
                                break;
                            }

                            // Go one step further in node path2
                            iter++;
                            // Find the next occurrence of the first element of
                            // node path 1.
                            iter =
                                std::find_if(iter, lastIndex + 1,
                                             [nodepath1](const uint32_t id) {
                                                 return id == (*nodepath1)[0];
                                             });
                        }
                    }
                }
            }
        }
    }

    // Iterate over all occurrences to determine whether they are non-redundant
    for (FMOccSFR& m : nonRedundantOcc) {
        if (!(m.isToBeRemoved() && !m.isToBeKept())) {
            // An occurrence is non-redundant if it is not the case that it is
            // marked as to be removed, but not as to be kept. An occurrence can
            // be to be removed as well as to be kept when its replacement does
            // not replace the occurrence it replaces. In this case, it is
            // non-redundant. An occurrence can also be marked as neither if it
            // has no other occurrence to be compared with. In this case, it is
            // also non-redundant.
            nonRedundantOccFinal.emplace_back(m);
        }
    }
}

template <class positionClass>
std::vector<FMOccSFR> FMIndexDBG<positionClass>::filterStrainFreeMatches(
    std::vector<FMOcc<FMPosSFR>>& occ, const int& maxED) {

    // Store whether linear or complete filtering must be applied.
    bool linear = !filteringOptionComplete;

    // Increase the total number of reported matches
    positionsInPostProcessingCounter += occ.size();

    for (auto& o : occ) {
        if (o.getPosition().getNodePath().empty()) {
            vector<uint32_t> nodePath;
            FMPosSFR pos = o.getPosition();
            findNodeUnderK(pos, 0, nodePath);
            pos.setNodePath(nodePath);
            o.setPosition(pos);
        }
    }

    // Sort the occurrences
    if (linear) {
        // For linear sorting: sort based on node path immediately, not taking
        // into account the node path length
        sort(occ.begin(), occ.end());
    } else {
        // For complete sorting: sort based on node path length first, then on
        // the node path itself
        sort(occ.begin(), occ.end(), compareComplete);
    }
    // Remove duplicates
    occ.erase(unique(occ.begin(), occ.end()), occ.end());

    if (occ.size() == 0) {
        // If there are no occurrences, return an empty list
        return {};
    } else if (occ.size() == 1) {
        // If there is only one occurrence, no filtering is needed. We transform
        // the occurrence type and return a vector with one element.
        vector<FMOccSFR> result;
        result.emplace_back(occ[0]);
        return result;
    }

    // Create the vector that will contain the occurrences without redundancy
    // for the same node path.
    std::vector<FMOccSFR> nonRedundantOcc;
    // Reserve space for the occurrences without redundancy for the same node
    // path.
    nonRedundantOcc.reserve(1000 * maxED);

    // Initialize the variables that will contain the minimum and maximum node
    // path lengths of all occurrences for this read
    size_t minLen, maxLen;
    // Filter all occurrences that have the same node path between themselves.
    // Add them to nonRedundantOcc. minLen and maxLen are also updated.
    // Additionally, the occurrences are transformed from FMOcc<FMPosSFR>
    // objects to FMOccSFR objects, which is necessary later on
    filterWithinSameNodePath(occ, nonRedundantOcc, maxED, minLen, maxLen);

    // Create the vector that will contain the final non-redundant occurrences.
    std::vector<FMOccSFR> nonRedundantOccFinal;
    // Reserve space for the non-redundant occurrences
    nonRedundantOccFinal.reserve(1000 * maxED);

    // Initiate the actual filtering process.
    if (linear) {
        filterDifferentNodePathsLinear(nonRedundantOcc, nonRedundantOccFinal,
                                       maxED, minLen, maxLen);
    } else {
        filterDifferentNodePathsComplete(nonRedundantOcc, nonRedundantOccFinal,
                                         maxED);
    }

    return nonRedundantOccFinal;
}

// ----------------------------------------------------------------------------
// HELPER ROUTINES FOR VISUALIZATION
// ----------------------------------------------------------------------------

template <class positionClass>
void FMIndexDBG<positionClass>::initializeFilesForVisualization(
    std::string filename, bool multipleSubgraphs, std::ofstream& edgefile) {
    // Initialize the file
    if (!multipleSubgraphs) {
        // Create empty file if this was not yet done
        edgefile.open(filename + "_SubgraphEdges.tsv");
    } else {
        // Append to existing file
        edgefile.open(filename + "_SubgraphEdges.tsv", std::ios_base::app);
    }

    if (!multipleSubgraphs) {
        // Write out the headers if this was not yet done
        edgefile
            << "EdgeKey\tSource\tOmegaShort\tOmegaFull\tPartOfPath\tTarget"
               "\tOmegaShort\tOmegaFull\tPartOfPath\tColor\tEdgeMultiplicity\n";
    }
}

template <class positionClass>
void FMIndexDBG<positionClass>::visitNode(
    uint32_t id, uint32_t depth, std::vector<uint32_t>& visited_nodes,
    std::queue<std::pair<int, int>>& node_queue) {
    // Push the node onto the queue
    node_queue.push(make_pair(id, depth));
    // Mark the node as visited
    G[id].visited = true;
    // Add the node to the visited nodes
    visited_nodes.emplace_back(id);
}

template <class positionClass>
void FMIndexDBG<positionClass>::fillInVisualizationNode(
    std::vector<VisualizationNode*>& visualizedNodes, uint32_t& id, Node& node,
    std::vector<uint32_t>& path) {

    // Check if this node is part of the original node path
    bool part_of_path = std::count(path.begin(), path.end(), id) > 0;
    // Get the index of an occurrence of this node in the text
    length_t indexInText = findSA(node.left_kmer_forward);
    // Get the string corresponding to the node
    string omega =
        text.decodeSubstring(sigma, indexInText, indexInText + node.len);
    string omega_short;

    // Get the short version of the string corresponding to the node
    if (omega.length() <= 15) {
        // The string corresponding to the node is short, so it is
        // visualized entirely
        omega_short = omega;
    } else {
        // The string corresponding to the node is long, so it is not
        // visualized entirely
        std::stringstream ss;
        ss << omega.substr(0, 5) << "..(+" << omega.length() - 10 << ").."
           << omega.substr(omega.length() - 5, 5);
        omega_short = ss.str();
    }

    // set the node pointer in the visualized nodes
    visualizedNodes[id] =
        new VisualizationNode(node, part_of_path, omega, omega_short);
}

template <class positionClass>
void FMIndexDBG<positionClass>::visualizeSubgraphIntermediary(
    std::vector<uint32_t>& path, std::string subgraph_id,
    std::vector<uint32_t>& visited_nodes,
    std::queue<std::pair<int, int>>& node_queue,
    std::map<std::string, std::map<length_t, length_t>>& edges,
    size_t& edgecounter, std::vector<VisualizationNode*>& visualizedNodes,
    bool separateEdges, std::set<uint32_t>& subgraphNodes) {

    // Take the next node ID from the queue along with its neighborhood depth
    std::pair<int, int> p = node_queue.front();
    node_queue.pop();
    // Find the corresponding node
    uint32_t id = p.first;
    Node node = G[id];
    subgraphNodes.insert(id);
    uint32_t current_depth = p.second;

    // Create the visualization node with extra attributes such as omega, if
    // necessary
    if (!visualizedNodes[id]) {
        fillInVisualizationNode(visualizedNodes, id, node, path);
    }

    // Store the visualization node
    VisualizationNode visNode = *visualizedNodes[id];

    // Iterate over all incoming edges and the corresponding predecessors
    for (size_t i = 0; i < node.multiplicity; i++) {
        uint32_t offset = i;
        uint32_t id_predecessor = jumpToPredecessorThroughEdge(id, offset);
        // Check that the corresponding predecessor is not an end node
        if (id_predecessor < numberOfGraphNodes - numberOfStrains) {
            // Store the predecessor node
            Node& nodePredecessor = G[id_predecessor];

            // If the predecessor has been visited, the edge needs to be
            // reported If the current depth is bigger than 0, the predecessor
            // needs to be visited as well
            if (nodePredecessor.visited || current_depth != 0) {

                // Create the visualization node with extra attributes such as
                // omega, if necessary
                if (!visualizedNodes[id_predecessor]) {
                    fillInVisualizationNode(visualizedNodes, id_predecessor,
                                            nodePredecessor, path);
                }

                // Store the visualization node
                VisualizationNode visNodePredecessor =
                    *visualizedNodes[id_predecessor];

                // Report the edge between the current node and its predecessor,
                // along with the node attributes to a temporary buffer
                std::stringstream buffer;
                buffer << subgraph_id << "_Edge" << edgecounter << "\t"
                       << subgraph_id << id_predecessor << "\t"
                       << id_predecessor << ":"
                       << visNodePredecessor.omega_short << "\t"
                       << visNodePredecessor.omega << "\t"
                       << visNodePredecessor.part_of_path << "\t" << subgraph_id
                       << id << "\t" << id << ":" << visNode.omega_short << "\t"
                       << visNode.omega << "\t" << visNode.part_of_path;
                if (separateEdges) {
                    // Each edge is reported separately with its corresponding
                    // strain ID
                    buffer << "\t"
                           << findStrain(findSA(i + node.left_kmer_forward));
                    edges[buffer.str()][0]++;
                } else {
                    // All edges between the same nodes are bundled into one.
                    // Strain IDs and multiplicities are stored as an attribute
                    edges[buffer.str()]
                         [findStrain(findSA(i + node.left_kmer_forward))]++;
                }
                // Make sure the predecessor is visited if it is within the
                // neighborhood depth
                if (!nodePredecessor.visited) {
                    visitNode(id_predecessor, current_depth - 1, visited_nodes,
                              node_queue);
                }
            }
        }
    }
    // Check that the current node is not an end note and that its successors
    // still belong to the neighborhood
    if (id < numberOfGraphNodes - numberOfStrains && current_depth > 0) {
        // Iterate over all possible characters
        for (size_t i = 0; i < sigma.size(); i++) {
            uint32_t id_successor;
            // Try to jump to a new character by appending a character to the
            // substring of the current node.
            if (jumpToSuccessorWithChar(id, id_successor, i)) {
                // A successor is present
                // Make sure the successor is visited
                if (!G[id_successor].visited) {
                    visitNode(id_successor, current_depth - 1, visited_nodes,
                              node_queue);
                }
                // If the character '%' was added, we must check if there are
                // multiple such end nodes
                if (i == 1) {
                    size_t offset = 1;
                    // Keep increasing the offset until no more successors are
                    // found
                    while (jumpToSuccessorWithChar(id, id_successor, i,
                                                   offset++)) {
                        // Make sure the successor is visited
                        if (!G[id_successor].visited) {
                            visitNode(id_successor, current_depth - 1,
                                      visited_nodes, node_queue);
                        }
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// ROUTINES FOR VISUALIZATION
// ----------------------------------------------------------------------------

template <class positionClass>
std::set<uint32_t> FMIndexDBG<positionClass>::visualizeSubgraph(
    std::vector<uint32_t>& path, uint32_t depth, std::string filename,
    bool separateEdges, bool multipleSubgraphs, string subgraph_id) {
    std::cout << "Constructing subgraph..." << std::endl;

    // Initialize the output file
    std::ofstream edgefile;
    initializeFilesForVisualization(filename, multipleSubgraphs, edgefile);

    // Save all nodes in the subgraph
    std::set<uint32_t> subgraphNodes;

    // Initialize a std::vector that keeps track of all visited nodes
    std::vector<uint32_t> visited_nodes{};
    // Create a node queue that keeps track of all nodes that still need to be
    // visited
    std::queue<std::pair<int, int>> node_queue;

    // Add each unique node in the node path to the node queue and visited nodes
    // once
    for (uint32_t id : path) {
        if (!G[id].visited) {
            visitNode(id, depth, visited_nodes, node_queue);
        }
    }

    // Initialize the edge counter that serves to create edge IDs
    size_t edgecounter = 0;

    // Create an empty array that will eventually contain all nodes visualized
    std::vector<VisualizationNode*> visualizedNodes(G.size(), nullptr);

    // Create an empty map that will eventually contain all edges visualized
    std::map<std::string, std::map<length_t, length_t>> edges;

    // Iterate over all nodes in the neighborhood of the node path
    while (!node_queue.empty()) {
        visualizeSubgraphIntermediary(
            path, subgraph_id, visited_nodes, node_queue, edges, edgecounter,
            visualizedNodes, separateEdges, subgraphNodes);
    }

    // Write out the edges to the actual output file
    for (const auto& edgePair : edges) {
        // Write out all attributes except for the strain ID and multiplicity
        // data
        edgefile << edgePair.first << "\t";
        // Initialize bool necessary for bundled edges
        bool firstEntry = true;
        // Initialize the total edge multiplicity, also necessary for bundled
        // edges
        length_t edgeMultiplicity = 0;
        // Iterate over all <ID,multiplicity> pairs
        for (const auto& subPair : edgePair.second) {
            if (separateEdges) {
                // No bundled edges, so there is only one edge pair present
                edgefile << subPair.second;
            } else {
                // Bundled edges
                if (firstEntry) {
                    firstEntry = false;
                } else {
                    // Separation character between <ID,multiplicity> pairs
                    edgefile << ",";
                }
                // Write out the strain ID with its multiplicity as color
                edgefile << subPair.first << ":" << subPair.second;
                // Keep track of the total multiplicity
                edgeMultiplicity += subPair.second;
            }
        }
        if (!separateEdges) {
            // For bundled edges: write out the total multiplicity
            edgefile << "\t" << edgeMultiplicity;
        }
        edgefile << "\n";
    }

    // Close output file
    edgefile.close();

    // Memory cleanup
    for (const auto& visNode : visualizedNodes) {
        if (visNode) {
            delete visNode;
        }
    }
    visualizedNodes.clear();

    // Reset the visited attribute of all visited nodes to enable the correct
    // construction of later subgraphs
    for (size_t i = 0; i < visited_nodes.size(); i++) {
        G[visited_nodes[i]].visited = false;
    }

    return subgraphNodes;
}

// For strain-fixed matching:
template <class positionClass>
void FMIndexDBG<positionClass>::visualizeSubgraphs(
    const std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>&
        paths,
    uint32_t depth, std::string filename, bool separateEdges) {
    std::ofstream edgefile;
    std::ofstream overviewfile;
    getText();
    // Create the output files
    edgefile.open(filename + "_SubgraphEdges.tsv");
    overviewfile.open(filename + "_SubgraphOverview.tsv");
    // Initialize the headers of the output files
    edgefile << "EdgeKey\tSource\tOmegaShort\tOmegaFull\tPartOfPath\tTarget\tOm"
                "egaShort\tOmegaFull\tPartOfPath\tColor\tEdgeMultiplicity\n";
    overviewfile << "SubgraphID\tPath\tDistanceFromLeftEnd\tStrain\tPosition\tL"
                    "ength\tED\n";

    edgefile.close();

    // Iterate over all paths corresponding to the matches of the pattern
    uint32_t counter = 0;
    for (const auto& path : paths) {
        // Visualize the subgraph corresponding to this path
        std::vector<uint32_t> nodepath = path.first;
        visualizeSubgraph(nodepath, depth, filename, separateEdges, true,
                          "Subgraph" + to_string(counter) + "_");
        // Report all matches in the reference text that correspond to this path
        for (length_t i = 0; i < path.second.size(); i++) {
            // Determine the separation character based on the occurrence length
            char separationchar =
                path.second[i].getRange().width() < k_DBG ? '/' : ',';
            // Report the path
            overviewfile << counter << "\t" << nodepath[0];
            for (length_t i = 1; i < nodepath.size(); i++) {
                overviewfile << separationchar << nodepath[i];
            }
            // Report information on the match in the reference
            overviewfile << "\t" << path.second[i].getDistanceFromLeftEnd()
                         << "\t" << path.second[i].getStrain() << "\t"
                         << path.second[i].getRange().getBegin() << "\t"
                         << path.second[i].getRange().width() << "\t"
                         << path.second[i].getDistance() << "\n";
        }
        counter++;
    }
    overviewfile.close();
}

// For strain-free matching:
template <class positionClass>
void FMIndexDBG<positionClass>::visualizeSubgraphs(
    const std::vector<FMOccSFR>& paths, uint32_t depth, std::string filename,
    bool separateEdges) {
    std::ofstream nodefile;
    std::ofstream edgefile;
    std::ofstream overviewfile;
    getText();
    // Create the output files
    edgefile.open(filename + "_SubgraphEdges.tsv");
    overviewfile.open(filename + "_SubgraphOverview.tsv");
    // Initialize the headers of the output files
    edgefile << "EdgeKey\tSource\tOmegaShort\tOmegaFull\tPartOfPath\tTarget\tOm"
                "egaShort\tOmegaFull\tPartOfPath\tColor\tEdgeMultiplicity\n";
    overviewfile << "SubgraphID\tPath\tDistanceFromLeftEnd\tLength\tED\n";

    edgefile.close();

    // Iterate over all paths corresponding to the matches of the pattern
    uint32_t counter = 0;
    for (const auto& path : paths) {
        // Visualize the subgraph corresponding to this path
        std::vector<uint32_t> nodepath = path.getPosition().getNodePath();
        visualizeSubgraph(nodepath, depth, filename, separateEdges, true,
                          "Subgraph" + to_string(counter) + "_");
        // Report the match in the reference text that corresponds to this path
        overviewfile << counter << "\t" << nodepath[0];
        for (length_t i = 1; i < nodepath.size(); i++) {
            overviewfile << "," << nodepath[i];
        }
        overviewfile << "\t" << path.getPosition().getDistanceFromLeftEnd()
                     << "\t" << path.getPosition().getTrueDepth() << "\t"
                     << path.getDistance() << "\n";

        counter++;
    }
    overviewfile.close();
}

template class FMIndexDBG<FMPos>;
template class FMIndexDBG<FMPosSFR>;