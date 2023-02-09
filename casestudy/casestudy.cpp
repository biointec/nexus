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

#include "../src/benchmarking.h"
#include "../src/searchstrategy.h"
#include "../src/strainfreemapper.h"
#include <algorithm>
#include <set>

vector<string> schemes = {"kuch1",  "kuch2", "kianfar", "manbest",
                          "pigeon", "01*0",  "custom",  "naive"};

struct neighborNode {
    // The name of the corresponding RRDR mutation
    string mutationName;
    // The identifier of the node
    uint32_t id;
    // The number of strains that pass through this node
    size_t numberOfMutationStrains;
    // The position of the possible compensatory mutation with respect to the
    // reference strain
    length_t position;
    // The length of the node
    uint32_t len;

    neighborNode(string mutationName, uint32_t id,
                 size_t numberOfMutationStrains, uint32_t len)
        : mutationName(mutationName), id(id),
          numberOfMutationStrains(numberOfMutationStrains), position(0),
          len(len) {
    }

    /**
     * @brief Writes the node to a file
     *
     * @param ofs output stream
     * @return std::ofstream& - output stream
     */
    std::ofstream& write(std::ofstream& ofs) {
        ofs << mutationName << "\t" << id << "\t" << numberOfMutationStrains
            << "\t" << position << "\t" << len << std::endl;

        return ofs;
    }
};

void handleMutation(const size_t& k, const FMIndexDBG<FMPos>& bwt,
                    const SearchStrategyDBG<FMIndexDBG<FMPos>, FMPos>* strategy,
                    std::set<uint32_t> neighboringNodesOnly,
                    const Range& rpoBRange, const Range& rpoCRange,
                    const Range& rpoARange, const Range& RRDR_range,
                    const string& mutationName, const string& mutationSequence,
                    const uint32_t& offsetOld, std::ofstream& outputFile) {
    // Get the node vector corresponding to the implicit graph
    std::vector<Node> G = bwt.getGraph();

    // Create a vector that will contain all strains that are mutated with the
    // current RRDR mutation
    std::set<uint32_t> mutationStrains;

    std::cout << "Investigating mutation " << mutationName << "..."
              << std::endl;

    // Sanity check
    if (mutationSequence.length() != k) {
        throw runtime_error("Mutation read should have length k.");
    }

    // Get the node that corresponds to the RRDR mutation.
    auto results = strategy->matchApproxSFI(mutationSequence, 0);
    auto it = results.begin();
    uint32_t nodeID = it->first[0];

    // Remove the node containing the mutation from the nodes that are to be
    // investigated. The RRDR mutation itself is not a compensatory mutation
    neighboringNodesOnly.erase(nodeID);

    // Insert all strains that contain the RRDR mutation in the designated
    // vector
    for (auto occ : it->second) {
        int strain = occ.getStrain();
        mutationStrains.insert(strain);
    }

    // Initialize a table that will contain compensatory mutations
    std::vector<neighborNode> table;

    // Iterate over all nodes in the remaining neighborhood (i.e.,
    // non-reference, not the RRDR mutation itself) and investigate their
    // passing strains
    for (uint32_t neighborID : neighboringNodesOnly) {
        // Extract the neighboring node
        Node neighbor = G[neighborID];
        // Extract its multiplicity
        uint32_t mult = neighbor.multiplicity;
        // Initialize a vector to store the strains that pass through this
        // neighboring node
        std::set<uint32_t> neighborStrains;
        // Iterate over all passing edges to find their strain IDs
        for (size_t i = 0; i < mult; i++) {
            // Add the strain to the designated set
            neighborStrains.insert(
                bwt.findStrain(bwt.findSA(neighbor.left_kmer_forward + i)));
        }
        // Find the strains that carry the RRDR mutation of interest AND that
        // pass through the neighboring node
        std::set<uint32_t> intersect;
        set_intersection(neighborStrains.begin(), neighborStrains.end(),
                         mutationStrains.begin(), mutationStrains.end(),
                         std::inserter(intersect, intersect.begin()));
        // Store the number of strains that carry the RRDR mutation of interest
        // AND that pass through the neighboring node
        size_t absNrOfMutated = intersect.size();
        // Check if the neighboring node meets our conditions: it must be
        // observed in at least two strains and it must carry ONLY strains that
        // also carry the RRDR mutation of interest
        if (absNrOfMutated > 1 && absNrOfMutated == neighborStrains.size()) {
            // If the node passes the conditions, create an entry and store it
            // in the table
            neighborNode nn(mutationName, neighborID, absNrOfMutated,
                            neighbor.len);
            table.push_back(nn);
        }
    }

    // Iterate over all entries in the table in order to try and find a
    // corresponding coordinate with respect to the reference strain (i.e., the
    // first strain in the concatenation)
    for (neighborNode& n : table) {
        // Initialize a vector that will contain the set of possible positions
        // corresponding to the presumed mutation in the node
        vector<length_t> possiblePositions;
        // Store nodes that have been visited in the process along with a
        // boolean indication whether they are a reference node (i.e., the
        // reference strain passes through them)
        map<uint32_t, bool> visitedNodesReference;
        // Iterate over all edges that pass through the possibly compensatory
        // node, to follow them towards a reference node
        for (size_t offset = 0; offset < G[n.id].multiplicity; offset++) {
            // This variable will contain identifiers for the predecessor node.
            // Initialize it with the neighboring node itself for now
            uint32_t predID = n.id;
            // Store the offset corresponding to the edge we are following in
            // this iteration
            uint32_t offsetToUpdate = offset;
            // Initialize a vector that will contain the set of possible
            // positions corresponding to the set of edges we are following here
            vector<length_t> positions;
            // We assume that k-1 character precede the actual mutated
            // nucleotide in this node (i.e., the k-1 characters that are shared
            // with the previous node)
            uint32_t extraDistance = k - 1;
            // We go back through the predecessors. We take at most 10 steps,
            // after which we decide that no position was found here.
            for (int i = 0; i < 10; i++) {
                // Find the predecessor ID by following the edges/strain in
                // upstream direction
                predID =
                    bwt.jumpToPredecessorThroughEdge(predID, offsetToUpdate);
                // Update the distance of the start of this node to the presumed
                // mutation
                extraDistance += G[predID].len - k + 1;
                // Check if this node was visited before
                if (visitedNodesReference.find(predID) ==
                    visitedNodesReference.end()) {
                    // The node was NOT visited before
                    // Iterate over all strains in the newly found node
                    for (size_t offset2 = 0; offset2 < G[predID].multiplicity;
                         offset2++) {
                        // Compute a position for the presumed compensatory
                        // mutation, by finding the position of this predecessor
                        // and adding the distance that we have been keeping
                        // along the way
                        length_t position =
                            bwt.findSA(G[predID].left_kmer_forward + offset2) +
                            extraDistance;
                        // Check two things: the position must be smaller than
                        // the text length, otherwise it is in valid. The
                        // position must also be part of the reference strain
                        // (i.e., the first one in the concatenation), because
                        // we want to assign a coordinate with respect to this
                        // reference strain.
                        if (position < bwt.getTextLength() &&
                            bwt.findStrain(position) == 0) {
                            // Conditions were met: add the position to the list
                            // of possible positions
                            positions.push_back(position);
                        }
                    }

                    if (positions.size() == 1) {
                        // Only one position was found, so the reference node
                        // tells us unambiguously where in the reference we are.
                        // We indicate it as usable in the vector
                        visitedNodesReference[predID] = true;
                        // Stop following this path, a position was found
                        break;
                    }
                    // The node is not a reference node OR multiple positions
                    // were found (ambiguity). So we clear the positions
                    positions.clear();
                    // Indicate that this not cannot be used to find a position
                    visitedNodesReference[predID] = false;
                } else if (visitedNodesReference[predID]) {
                    // The node was visited before and can be used to
                    // unambiguously locate. It has already been used, so no
                    // more work needs to be done here on this path.
                    break;
                }
            }

            // Check if a any positions were found for this edge in the possibly
            // compensatory node
            if (!positions.empty()) {
                // Iterate over all positions
                for (auto position : positions) {
                    // Check that this position is related to one of the three
                    // genes of interest. Otherwise, it was found in a node that
                    // coincidentally also lies in the neighborhood, but this is
                    // not of interest.
                    if (rpoARange.contains(position) ||
                        rpoBRange.contains(position) ||
                        rpoCRange.contains(position)) {
                        // Position is inside one of the three genes, so it
                        // still qualifies as a possible position.
                        possiblePositions.push_back(position);
                    }
                }
            }
        }
        // Initialize the final position
        length_t finalPosition = 0;
        // Check if there are any possibilities
        if (!possiblePositions.empty()) {
            // Choose the maximum element, which is assumed to be found through
            // the shortest path from mutation to reference
            finalPosition = *max_element(possiblePositions.begin(),
                                         possiblePositions.end());
        }
        // store the possible compensatory nodes final position
        n.position = finalPosition;
    }

    // Iterate again over all possible compensatory nodes
    for (neighborNode n : table) {
        // Check that their position was initialized and is not within the RRDR
        // region (which is not where we are looking)
        if (n.position && !RRDR_range.contains(n.position)) {
            // Write the candidate putative compensatory mutation to the output
            // file
            n.write(outputFile);
        }
    }
}

int main(int argc, char* argv[]) {

    // De Bruijn parameter is 19, which was manually decided to be the best
    // option
    size_t k = 19;

    // The filename for the pan-genome. This program must be run in the folder
    // that contains the file.
    std::string filename = "MTuberculosisPanGenome";

    // Reconstruct the index
    FMIndexDBG<FMPos> bwt(filename, 16, 128, k, false);

    // Construct a search strategy to search the graph
    SearchStrategyDBG<FMIndexDBG<FMPos>, FMPos>* strategy;
    strategy = new KucherovKplus1DBG<FMIndexDBG<FMPos>, FMPos>(bwt, DYNAMIC,
                                                               EDITOPTIMIZED);

    // Store the sequence corresponding to the rpoB gene. This is obtained from
    // here:
    // https://www.ncbi.nlm.nih.gov/nuccore/CP003248.2?report=fasta&from=759810&to=763328
    string rpoBSequence =
        "TTGGCAGATTCCCGCCAGAGCAAAACAGCCGCTAGTCCTAGTCCGAGTCGCCCGCAAAGTTCCTCGAATA"
        "ACTCCGTACCCGGAGCGCCAAACCGGGTCTCCTTCGCTAAGCTGCGCGAACCACTTGAGGTTCCGGGACT"
        "CCTTGACGTCCAGACCGATTCGTTCGAGTGGCTGATCGGTTCGCCGCGCTGGCGCGAATCCGCCGCCGAG"
        "CGGGGTGATGTCAACCCAGTGGGTGGCCTGGAAGAGGTGCTCTACGAGCTGTCTCCGATCGAGGACTTCT"
        "CCGGGTCGATGTCGTTGTCGTTCTCTGACCCTCGTTTCGACGATGTCAAGGCACCCGTCGACGAGTGCAA"
        "AGACAAGGACATGACGTACGCGGCTCCACTGTTCGTCACCGCCGAGTTCATCAACAACAACACCGGTGAG"
        "ATCAAGAGTCAGACGGTGTTCATGGGTGACTTCCCGATGATGACCGAGAAGGGCACGTTCATCATCAACG"
        "GGACCGAGCGTGTGGTGGTCAGCCAGCTGGTGCGGTCGCCCGGGGTGTACTTCGACGAGACCATTGACAA"
        "GTCCACCGACAAGACGCTGCACAGCGTCAAGGTGATCCCGAGCCGCGGCGCGTGGCTCGAGTTTGACGTC"
        "GACAAGCGCGACACCGTCGGCGTGCGCATCGACCGCAAACGCCGGCAACCGGTCACCGTGCTGCTCAAGG"
        "CGCTGGGCTGGACCAGCGAGCAGATTGTCGAGCGGTTCGGGTTCTCCGAGATCATGCGATCGACGCTGGA"
        "GAAGGACAACACCGTCGGCACCGACGAGGCGCTGTTGGACATCTACCGCAAGCTGCGTCCGGGCGAGCCC"
        "CCGACCAAAGAGTCAGCGCAGACGCTGTTGGAAAACTTGTTCTTCAAGGAGAAGCGCTACGACCTGGCCC"
        "GCGTCGGTCGCTATAAGGTCAACAAGAAGCTCGGGCTGCATGTCGGCGAGCCCATCACGTCGTCGACGCT"
        "GACCGAAGAAGACGTCGTGGCCACCATCGAATATCTGGTCCGCTTGCACGAGGGTCAGACCACGATGACC"
        "GTTCCGGGCGGCGTCGAGGTGCCGGTGGAAACCGACGACATCGACCACTTCGGCAACCGCCGCCTGCGTA"
        "CGGTCGGCGAGCTGATCCAAAACCAGATCCGGGTCGGCATGTCGCGGATGGAGCGGGTGGTCCGGGAGCG"
        "GATGACCACCCAGGACGTGGAGGCGATCACACCGCAGACGTTGATCAACATCCGGCCGGTGGTCGCCGCG"
        "ATCAAGGAGTTCTTCGGCACCAGCCAGCTGAGCCAATTCATGGACCAGAACAACCCGCTGTCGGGGTTGA"
        "CCCACAAGCGCCGACTGTCGGCGCTGGGGCCCGGCGGTCTGTCACGTGAGCGTGCCGGGCTGGAGGTCCG"
        "CGACGTGCACCCGTCGCACTACGGCCGGATGTGCCCGATCGAAACCCCTGAGGGGCCCAACATCGGTCTG"
        "ATCGGCTCGCTGTCGGTGTACGCGCGGGTCAACCCGTTCGGGTTCATCGAAACGCCGTACCGCAAGGTGG"
        "TCGACGGCGTGGTTAGCGACGAGATCGTGTACCTGACCGCCGACGAGGAGGACCGCCACGTGGTGGCACA"
        "GGCCAATTCGCCGATCGATGCGGACGGTCGCTTCGTCGAGCCGCGCGTGCTGGTCCGCCGCAAGGCGGGC"
        "GAGGTGGAGTACGTGCCCTCGTCTGAGGTGGACTACATGGACGTCTCGCCCCGCCAGATGGTGTCGGTGG"
        "CCACCGCGATGATTCCCTTCCTGGAGCACGACGACGCCAACCGTGCCCTCATGGGGGCAAACATGCAGCG"
        "CCAGGCGGTGCCGCTGGTCCGTAGCGAGGCCCCGCTGGTGGGCACCGGGATGGAGCTGCGCGCGGCGATC"
        "GACGCCGGCGACGTCGTCGTCGCCGAAGAAAGCGGCGTCATCGAGGAGGTGTCGGCCGACTACATCACTG"
        "TGATGCACGACAACGGCACCCGGCGTACCTACCGGATGCGCAAGTTTGCCCGGTCCAACCACGGCACTTG"
        "CGCCAACCAGTGCCCCATCGTGGACGCGGGCGACCGAGTCGAGGCCGGTCAGGTGATCGCCGACGGTCCC"
        "TGTACTGACGACGGCGAGATGGCGCTGGGCAAGAACCTGCTGGTGGCCATCATGCCGTGGGAGGGCCACA"
        "ACTACGAGGACGCGATCATCCTGTCCAACCGCCTGGTCGAAGAGGACGTGCTCACCTCGATCCACATCGA"
        "GGAGCATGAGATCGATGCTCGCGACACCAAGCTGGGTGCGGAGGAGATCACCCGCGACATCCCGAACATC"
        "TCCGACGAGGTGCTCGCCGACCTGGATGAGCGGGGCATCGTGCGCATCGGTGCCGAGGTTCGCGACGGGG"
        "ACATCCTGGTCGGCAAGGTCACCCCGAAGGGTGAGACCGAGCTGACGCCGGAGGAGCGGCTGCTGCGTGC"
        "CATCTTCGGTGAGAAGGCCCGCGAGGTGCGCGACACTTCGCTGAAGGTGCCGCACGGCGAATCCGGCAAG"
        "GTGATCGGCATTCGGGTGTTTTCCCGCGAGGACGAGGACGAGTTGCCGGCCGGTGTCAACGAGCTGGTGC"
        "GTGTGTATGTGGCTCAGAAACGCAAGATCTCCGACGGTGACAAGCTGGCCGGCCGGCACGGCAACAAGGG"
        "CGTGATCGGCAAGATCCTGCCGGTTGAGGACATGCCGTTCCTTGCCGACGGCACCCCGGTGGACATTATT"
        "TTGAACACCCACGGCGTGCCGCGACGGATGAACATCGGCCAGATTTTGGAGACCCACCTGGGTTGGTGTG"
        "CCCACAGCGGCTGGAAGGTCGACGCCGCCAAGGGGGTTCCGGACTGGGCCGCCAGGCTGCCCGACGAACT"
        "GCTCGAGGCGCAGCCGAACGCCATTGTGTCGACGCCGGTGTTCGACGGCGCCCAGGAGGCCGAGCTGCAG"
        "GGCCTGTTGTCGTGCACGCTGCCCAACCGCGACGGTGACGTGCTGGTCGACGCCGACGGCAAGGCCATGC"
        "TCTTCGACGGGCGCAGCGGCGAGCCGTTCCCGTACCCGGTCACGGTTGGCTACATGTACATCATGAAGCT"
        "GCACCACCTGGTGGACGACAAGATCCACGCCCGCTCCACCGGGCCGTACTCGATGATCACCCAGCAGCCG"
        "CTGGGCGGTAAGGCGCAGTTCGGTGGCCAGCGGTTCGGGGAGATGGAGTGCTGGGCCATGCAGGCCTACG"
        "GTGCTGCCTACACCCTGCAGGAGCTGTTGACCATCAAGTCCGATGACACCGTCGGCCGCGTCAAGGTGTA"
        "CGAGGCGATCGTCAAGGGTGAGAACATCCCGGAGCCGGGCATCCCCGAGTCGTTCAAGGTGCTGCTCAAA"
        "GAACTGCAGTCGCTGTGCCTCAACGTCGAGGTGCTATCGAGTGACGGTGCGGCGATCGAACTGCGCGAAG"
        "GTGAGGACGAGGACCTGGAGCGGGCCGCGGCCAACCTGGGAATCAATCTGTCCCGCAACGAATCCGCAAG"
        "TGTCGAGGATCTTGCGTAA";
    // Store the corresponding coordinate range found at NCBI (gene RVBD_0667)
    // Zero-based indexing
    Range rpoBRange(759809, 763328);

    // Store the sequence corresponding to the rpoB gene. This is obtained from
    // here:
    // https://www.ncbi.nlm.nih.gov/nuccore/CP003248.2?report=fasta&from=763373&to=767323
    string rpoCSequence =
        "GTGCTCGACGTCAACTTCTTCGATGAACTCCGCATCGGTCTTGCTACCGCGGAGGACATCAGGCAATGGT"
        "CCTATGGCGAGGTCAAAAAGCCGGAGACGATCAACTACCGCACGCTTAAGCCGGAGAAGGACGGCCTGTT"
        "CTGCGAGAAGATCTTCGGGCCGACTCGCGACTGGGAATGCTACTGCGGCAAGTACAAGCGGGTGCGCTTC"
        "AAGGGCATCATCTGCGAGCGCTGCGGCGTCGAGGTGACCCGCGCCAAGGTGCGTCGTGAGCGGATGGGCC"
        "ACATCGAGCTTGCCGCGCCCGTCACCCACATCTGGTACTTCAAGGGTGTGCCCTCGCGGCTGGGGTATCT"
        "GCTGGACCTGGCCCCGAAGGACCTGGAGAAGATCATCTACTTCGCTGCCTACGTGATCACCTCGGTCGAC"
        "GAGGAGATGCGCCACAATGAGCTCTCCACGCTCGAGGCCGAAATGGCGGTGGAGCGCAAGGCCGTCGAAG"
        "ACCAGCGCGACGGCGAACTAGAGGCCCGGGCGCAAAAGCTGGAGGCCGACCTGGCCGAGCTGGAGGCCGA"
        "GGGCGCCAAGGCCGATGCGCGGCGCAAGGTTCGCGACGGCGGCGAGCGCGAGATGCGCCAGATCCGTGAC"
        "CGCGCGCAGCGTGAGCTGGACCGGTTGGAGGACATCTGGAGCACTTTCACCAAGCTGGCGCCCAAGCAGC"
        "TGATCGTCGACGAAAACCTCTACCGCGAACTCGTCGACCGCTACGGCGAGTACTTCACCGGTGCCATGGG"
        "CGCGGAGTCGATCCAGAAGCTGATCGAGAACTTCGACATCGACGCCGAAGCCGAGTCGCTGCGGGATGTC"
        "ATCCGAAACGGCAAGGGGCAGAAGAAGCTTCGCGCCCTCAAGCGGCTGAAGGTGGTTGCGGCGTTCCAAC"
        "AGTCGGGCAACTCGCCGATGGGCATGGTGCTCGACGCCGTCCCGGTGATCCCGCCGGAGCTGCGCCCGAT"
        "GGTGCAGCTCGACGGCGGCCGGTTCGCCACGTCCGACTTGAACGACCTGTACCGCAGGGTGATCAACCGC"
        "AACAACCGGCTGAAAAGGCTGATCGATCTGGGTGCGCCGGAAATCATCGTCAACAACGAGAAGCGGATGC"
        "TGCAGGAATCCGTGGACGCGCTGTTCGACAATGGCCGCCGCGGCCGGCCCGTCACCGGGCCGGGCAACCG"
        "TCCGCTCAAGTCGCTTTCCGATCTGCTCAAGGGCAAGCAGGGCCGGTTCCGGCAGAACCTGCTCGGCAAG"
        "CGTGTCGACTACTCGGGCCGGTCGGTCATCGTGGTCGGCCCGCAGCTCAAGCTGCACCAGTGCGGTCTGC"
        "CCAAGCTGATGGCGCTGGAGCTGTTCAAGCCGTTCGTGATGAAGCGGCTGGTGGACCTCAACCATGCGCA"
        "GAACATCAAGAGCGCCAAGCGCATGGTGGAGCGCCAGCGCCCCCAAGTGTGGGATGTGCTCGAAGAGGTC"
        "ATCGCCGAGCACCCGGTGTTGCTGAACCGCGCACCCACCCTGCACCGGTTGGGTATCCAGGCCTTCGAGC"
        "CAATGCTGGTGGAAGGCAAGGCCATTCAGCTGCACCCGTTGGTGTGTGAGGCGTTCAATGCCGACTTCGA"
        "CGGTGACCAGATGGCCGTGCACCTGCCTTTGAGCGCCGAAGCGCAGGCCGAGGCTCGCATTTTGATGTTG"
        "TCCTCCAACAACATCCTGTCGCCGGCATCTGGGCGTCCGTTGGCCATGCCGCGGCTGGACATGGTGACCG"
        "GGCTGTACTACCTGACCACCGAGGTCCCCGGGGACACCGGCGAATACCAGCCGGCCAGCGGGGATCACCC"
        "GGAGACTGGTGTCTACTCTTCGCCGGCCGAAGCGATCATGGCGGCCGACCGCGGTGTCTTGAGCGTGCGG"
        "GCCAAGATCAAGGTGCGGCTGACCCAGCTGCGGCCGCCGGTCGAGATCGAGGCCGAGCTATTCGGCCACA"
        "GCGGCTGGCAGCCGGGCGATGCGTGGATGGCCGAGACCACGCTGGGCCGGGTGATGTTCAACGAGCTGCT"
        "GCCGCTGGGTTATCCGTTCGTCAACAAGCAGATGCACAAGAAGGTGCAGGCCGCCATCATCAACGACCTG"
        "GCCGAGCGTTACCCGATGATCGTGGTCGCCCAGACCGTCGACAAGCTCAAGGACGCCGGCTTCTACTGGG"
        "CCACCCGCAGCGGCGTGACGGTGTCGATGGCCGACGTGCTGGTGCCGCCGCGCAAGAAGGAGATCCTCGA"
        "CCACTACGAGGAGCGCGCGGACAAGGTCGAAAAGCAGTTCCAGCGTGGCGCTTTGAACCACGACGAGCGC"
        "AACGAGGCGCTGGTGGAGATTTGGAAGGAAGCCACCGACGAGGTCGGTCAGGCGTTGCGGGAGCACTACC"
        "CCGACGACAACCCGATCATCACCATCGTCGACTCCGGCGCCACCGGCAACTTCACCCAGACTCGAACGCT"
        "GGCCGGTATGAAGGGCCTGGTGACCAACCCGAAGGGTGAGTTCATCCCGCGTCCGGTCAAGTCCTCCTTC"
        "CGTGAGGGCCTGACCGTGCTGGAGTACTTCATCAACACCCACGGCGCTCGAAAGGGCTTGGCGGACACCG"
        "CGTTGCGCACCGCCGACTCCGGCTACCTGACCCGACGTCTGGTGGACGTGTCCCAGGACGTGATCGTGCG"
        "CGAGCACGACTGCCAGACCGAGCGCGGCATCGTCGTCGAGCTGGCCGAGCGTGCACCCGACGGCACGCTG"
        "ATCCGCGACCCGTACATCGAAACCTCGGCCTACGCGCGGACCCTGGGCACCGACGCGGTCGACGAGGCCG"
        "GCAACGTCATCGTCGAGCGTGGTCAAGACCTGGGCGATCCGGAGATTGACGCTCTGTTGGCTGCTGGTAT"
        "TACCCAGGTCAAGGTGCGTTCGGTGCTGACGTGTGCCACCAGCACCGGCGTGTGCGCGACCTGCTACGGG"
        "CGTTCCATGGCCACCGGCAAGCTGGTCGACATCGGTGAAGCCGTCGGCATCGTGGCCGCCCAGTCCATCG"
        "GCGAACCCGGCACCCAGCTGACCATGCGCACCTTCCACCAGGGTGGCGTCGGTGAGGACATCACCGGTGG"
        "TCTGCCCCGGGTGCAGGAGCTGTTCGAGGCCCGGGTACCGCGTGGCAAGGCGCCGATCGCCGACGTCACC"
        "GGCCGGGTTCGGCTCGAGGACGGCGAGCGGTTCTACAAGATCACCATCGTTCCTGACGACGGCGGTGAGG"
        "AAGTGGTCTACGACAAGATCTCCAAGCGGCAGCGGCTGCGGGTGTTCAAGCACGAAGACGGTTCCGAACG"
        "GGTGCTCTCCGATGGCGACCACGTCGAGGTGGGCCAGCAGCTGATGGAAGGCTCGGCCGACCCGCATGAG"
        "GTGCTGCGGGTGCAGGGCCCCCGCGAGGTGCAGATACACCTGGTTCGCGAGGTCCAGGAGGTCTACCGCG"
        "CCCAAGGTGTGTCGATCCACGACAAGCACATCGAGGTGATCGTTCGCCAGATGCTGCGCCGGGTGACCAT"
        "CATCGACTCGGGCTCGACGGAGTTTTTGCCTGGCTCGCTGATCGACCGCGCGGAGTTCGAGGCAGAGAAC"
        "CGCCGAGTGGTGGCCGAGGGCGGTGAGCCCGCGGCCGGCCGTCCGGTGCTGATGGGCATCACGAAGGCGT"
        "CGCTGGCCACCGACTCGTGGCTGTCGGCGGCGTCGTTCCAGGAGACCACTCGCGTGCTGACCGATGCGGC"
        "GATCAACTGCCGCAGCGATAAGCTCAACGGTCTGAAGGAAAACGTGATCATCGGCAAGCTGATCCCGGCC"
        "GGTACCGGTATCAACCGCTACCGCAACATCGCGGTGCAGCCCACCGAGGAGGCCCGCGCTGCGGCGTACA"
        "CCATCCCGTCGTATGAGGATCAGTACTACAGCCCGGACTTCGGTGCGGCCACCGGTGCTGCCGTCCCGCT"
        "GGACGACTACGGCTACAGCGACTACCGCTAG";
    // Store the corresponding coordinate range found at NCBI (gene RVBD_0668)
    // Zero-based indexing
    Range rpoCRange(763372, 767323);

    // Store the sequence corresponding to the rpoB gene. This is obtained from
    // here:
    // https://www.ncbi.nlm.nih.gov/nuccore/CP003248.2?report=fasta&from=3877643&to=3878686
    string rpoASequence =
        "CTAAAGCTGTTCGGTTTCGGCGTAGTCCTGCTCGTCGTACGCGCCCTCGGTCGACCAGGTGCCGGTGGCG"
        "ACGTCGTAGCCCGCGACCTCCGAGGGGTCGAAGCTCGGCGGGCTGTCCTTGAGTGACAGGCCCAGCTGGT"
        "GCAGCTTGATCTTCACCTCGTCGATGGACTTCTGACCGAAGTTGCGGATGTCAAGCAGGTCGGATTCGGT"
        "GCGCGCCACCAGTTCGCCCACGGTGTGCACCCCCTCGCGCTTGAGGCAGTTGTAGGACCGCACCGTCAGA"
        "TCCAGGTCGTCGATCGGCAGGGCGAATGACGCAATGTGATCGGCCTCGGCCGGCGACGGCCCGATCTCGA"
        "TGCCTTCGGCCTCGACGTTGAGTTCCCGTGCCAGGCCGAACAACTCGACCAGCGTCTTGCCAGCCGACGC"
        "CAGCGCGTCGCGCGGGCTGATTGAATTCTTGGTCTCCACGTCCAGGATCAGCTTGTCGAAGTCGGTGCGC"
        "TGCTCGACCCGGGTGGCGTCCACCTTGTAGGTCACTTTGAGCACCGGTGAGTAGATGGAATCGACTGGAA"
        "TGCGCCCAATTTCGGCACCCGAAGCCCGGTTTTGCACCGCCGGGACATAGCCGCGGCCACGCTCGACGAC"
        "GAGCTCGACTTCCAGCTTGCCCTTATCGTTCAGCGTGGCGATGTGCATGCCGGGGTTGTGCACGGTGACG"
        "CCGGCCGGCGGCACGATGTCGCCGGCGGTAACCTCACCCGGACCCTGCTTGCGTAGGTACATGGTGACCG"
        "GCTCGTCCTCCTCCGAGGACACCACCAGGCTCTTGAGATTCAGGATGATCTCGGTGACATCTTCTTTGAC"
        "CCCGGGCACCGTGGTGAATTCGTGCAGTACACCATCGATGCGAATGCTGGTGACGGCCGCTCCGGGAATC"
        "GACGACAGCAGGGTGCGACGCAGCGAATTGCCCAGGGTGTAGCCGAATCCCGGCTCCAGCGGTTCGATCA"
        "CGAACTGGGATCGGTTGTCGGTGAGGACGTCCTCGGACAGGGTGGGGCGCTGTGAGATCAGCAT";
    // Store the corresponding coordinate range found at NCBI (gene RVBD_3457c)
    // Zero-based indexing
    Range rpoARange(3877642, 3878686);

    // Store the sequence corresponding to the RRDR region. This is obtained
    // from literature
    string RRDRSequence = "GGCACCAGCCAGCTGAGCCAATTCATGGACCAGAACAACCCGCTGTCGGGGT"
                          "TGACCCACAAGCGCCGACTGTCGGCGCTG";
    // Store the coordinate range for the RRDR region, obtained by matching the
    // sequence exactly to our reference
    // Zero-based indexing
    Range RRDRrange(761084, 761084 + 81);

    // Find the node paths corresponding to rpoB, rpoC, rpoA and RRDR
    auto resultsRpoB = strategy->matchApproxSFI(rpoBSequence, 0);
    auto resultsRpoC = strategy->matchApproxSFI(rpoCSequence, 0);
    auto resultsRpoA = strategy->matchApproxSFI(rpoASequence, 0);
    auto resultsRRDR = strategy->matchApproxSFI(RRDRSequence, 0);

    // Sanity check: each of the previous results should have exactly 1
    // occurrence in the graph, and 1 text occurrence in the reference strain
    // (ID 0)
    if (resultsRpoB.size() != 1 && resultsRpoC.size() != 1 &&
        resultsRpoA.size() != 1 && resultsRRDR.size() != 1) {
        throw runtime_error("There should be only 1 match in the graph.");
    }
    size_t referenceCounter = 0;
    for (auto occ : resultsRpoB.begin()->second) {
        if (occ.getStrain() == 0) {
            referenceCounter++;
        }
    }
    if (referenceCounter != 1) {
        throw runtime_error(
            "There should be only 1 match in the reference for rpoB.");
    }
    referenceCounter = 0;
    for (auto occ : resultsRpoC.begin()->second) {
        if (occ.getStrain() == 0) {
            referenceCounter++;
        }
    }
    if (referenceCounter != 1) {
        throw runtime_error(
            "There should be only 1 match in the reference for rpoC.");
    }
    referenceCounter = 0;
    for (auto occ : resultsRpoA.begin()->second) {
        if (occ.getStrain() == 0) {
            referenceCounter++;
        }
    }
    if (referenceCounter != 1) {
        throw runtime_error(
            "There should be only 1 match in the reference for rpoA.");
    }
    referenceCounter = 0;
    for (auto occ : resultsRRDR.begin()->second) {
        if (occ.getStrain() == 0) {
            referenceCounter++;
        }
    }
    if (referenceCounter != 1) {
        throw runtime_error(
            "There should be only 1 match in the reference for RRDR.");
    }

    // Create a vector containing the node paths corresponding to rpoA, rpoB and
    // rpoC. The neighborhood of these will be of interest for compensatory
    // mutations.
    std::vector<uint32_t> referenceNodePaths = resultsRpoB.begin()->first;
    referenceNodePaths.insert(referenceNodePaths.end(),
                              resultsRpoC.begin()->first.begin(),
                              resultsRpoC.begin()->first.end());
    referenceNodePaths.insert(referenceNodePaths.end(),
                              resultsRpoA.begin()->first.begin(),
                              resultsRpoA.begin()->first.end());

    // Define a maximum depth for nodes in the neighborhood which we will
    // consider to be possible compensatory mutations.
    int maxDepth = 5;
    std::cout << "Investigating all neighboring nodes at a maximum distance of "
              << maxDepth << "..." << std::endl;

    // Create the set of nodes in the neighborhood of our three genes of
    // interest
    std::set<uint32_t> neighborhoodNodes =
        bwt.visualizeSubgraph(referenceNodePaths, maxDepth, filename, false);

    // Create a set that will contain the neighboring nodes, without the nodes
    // of the reference path itself. Only nodes with a variation can possibly
    // compensatory mutations.
    std::set<uint32_t> neighboringNodesOnly;
    // Take the neighborhood minus the reference path
    std::set_difference(
        neighborhoodNodes.begin(), neighborhoodNodes.end(),
        referenceNodePaths.begin(), referenceNodePaths.end(),
        std::inserter(neighboringNodesOnly, neighboringNodesOnly.end()));

    std::cout << neighboringNodesOnly.size()
              << " possible compensatory nodes are being considered."
              << std::endl;

    // Store the three RRDR mutations from Table 12 in the manuscript for which
    // we try to find compensatory mutations For each mutation, we store the
    // mutation name, and the sequence that starts with the mutated nucleotide
    // and is padded with 18 additional characters from the reference genome in
    // the back. This way, we can assign a unique graph node to the mutation,
    // which is then the LAST graph node in which the mutation appears. In
    // reality, these three mutations each only appear in 1 node. If this had
    // not been the case, extra attention would have been necessary to avoid
    // false positive compensatory mutations etc.
    std::vector<std::tuple<string, string, uint32_t>> mutations = {
        std::make_tuple("D435G", "GCCAGAACAACCCGCTGTC", 0),
        std::make_tuple("S450L", "TGGCGCTGGGGCCCGGCGG", 0),
        std::make_tuple("L452P", "CGGGGCCCGGCGGTCTGTC", 0)};

    // Open an output file to which we will write the potential compensatory
    // mutations
    std::ofstream outputFile;
    outputFile.open(filename + "_Compensatory.tsv");

    // Write the header to the output file
    outputFile << "MutationName\tNodeID\tNumberOfSTrains\tpositionInReference\t"
                  "nodeLength\n";

    // Now execute compensatory mutation search process for each of the three
    // mutations of interest.
    for (auto p : mutations) {
        handleMutation(k, bwt, strategy, neighboringNodesOnly, rpoBRange,
                       rpoCRange, rpoARange, RRDRrange, std::get<0>(p),
                       std::get<1>(p), std::get<2>(p), outputFile);
    }

    outputFile.close();

    delete strategy;

    cout << "Bye...\n";
}
