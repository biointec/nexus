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

#include "benchmarking.h"
#include "strainfreemapper.h"

#include <chrono>  // timing
#include <numeric> // for summing over vector
#include <sstream> // used for splitting strings when reading csv files (outdated)

using namespace std;

vector<string> schemes = {"kuch1",  "kuch2", "kianfar", "manbest",
                          "pigeon", "01*0",  "custom",  "naive"};

string getFileExt(const string& s) {

    size_t i = s.rfind('.', s.length());
    if (i != string::npos) {
        return (s.substr(i + 1, s.length() - i));
    }

    return ("");
}

size_t getReads(vector<pair<string, string>>& reads, string& file,
                ifstream& ifile, size_t chunkSize, string& line,
                bool readWithN) {

    string read = "";
    string p = "";

    size_t chunkCounter = 0;
    bool readLine = false;

    // Read first line of chunk
    if (!line.empty() && (line[0] == '>' || line[0] == '@')) {
        p = (line.substr(1));
        readLine = true;
    }

    reads.reserve(2 * chunkSize);

    const auto& extension = getFileExt(file);

    bool fasta =
        (extension == "FASTA") || (extension == "fasta") || (extension == "fa");
    bool fastq = (extension == "fq") || (extension == "fastq");

    if (!ifile) {
        throw runtime_error("Cannot open file " + file);
    }
    if (!fasta && !fastq) {
        // this is a csv file
        if (extension != "csv") {
            throw runtime_error("extension " + extension +
                                " is not a valid extension for the readsfile");
        }
        string line;
        // get the first line we do not need this
        std::getline(ifile, line);

        while (std::getline(ifile, line)) {
            istringstream iss{line};

            vector<string> tokens;
            string token;

            while (getline(iss, token, ',')) {
                tokens.push_back(token);
            }
            string position = tokens[1];
            string read =
                tokens[2]; // ED + 2 column contains a read with ED compared
                           // to the read at position position with length
            string p = position;
            reads.push_back(make_pair(p, read));
        }
    } else if (fasta) {
        while (chunkCounter < chunkSize && std::getline(ifile, line)) {
            if (!line.empty() && line[0] == '>') {

                if (!read.empty()) {
                    // Ignore reads containing N
                    if (read.find('N') == std::string::npos) {
                        reads.push_back(make_pair(p, read));
                        reads.push_back(
                            make_pair(p, Nucleotide::getRevCompl(read)));
                        read.clear();
                        chunkCounter++;
                    } else {
                        readWithN = true;
                    }
                }

                p = (line.substr(1));

            } else {
                read += line;
            }
        }
        if (!read.empty()) {
            // Ignore reads containing N
            if (read.find('N') == std::string::npos) {
                reads.push_back(make_pair(p, read));
                reads.push_back(make_pair(p, Nucleotide::getRevCompl(read)));
                read.clear();
                chunkCounter++;
            } else {
                readWithN = true;
            }
        }
    } else {
        // fastQ
        while (chunkCounter < chunkSize && std::getline(ifile, line)) {
            if (!line.empty() && line[0] == '@') {
                if (!read.empty()) {
                    if (read.find('N') == std::string::npos) {
                        reads.push_back(make_pair(p, read));
                        reads.push_back(
                            make_pair(p, Nucleotide::getRevCompl(read)));
                        read.clear();
                        chunkCounter++;
                    } else {
                        readWithN = true;
                    }
                }
                p = (line.substr(1));
                readLine = true;
            } else if (readLine) {
                read = line;
                readLine = false;
            }
        }
        if (!read.empty()) {
            if (read.find('N') == std::string::npos) {
                reads.push_back(make_pair(p, read));
                reads.push_back(make_pair(p, Nucleotide::getRevCompl(read)));
                read.clear();
                chunkCounter++;
            } else {
                readWithN = true;
            }
        }
    }
    return chunkCounter;
}

double
avgVec(vector<length_t> const& v) // note: the average must not be an integer
{
    return v.empty() ? 0.0 : accumulate(v.begin(), v.end(), 0.0) / v.size();
    ;
}

length_t sum(vector<length_t> const& v) {
    return accumulate(v.begin(), v.end(), 0.0);
}

void writeToOutputSFI(
    const string& file,
    const vector<std::map<std::vector<uint32_t>,
                          std::vector<TextOccurrenceSFI>>>& mPerRead,
    const vector<pair<string, string>>& reads, bool& firstChunk, ofstream& f2) {

    if (firstChunk) {
        // Write header
        firstChunk = false;

        // cout << "Writing to output file " << file << " ..." << endl;

        f2 << "Identifier\tSubgraphID\tPath\tDistanceFromLeftEnd\tStrain\t"
              "Position\tLength\tED\treverseComplement\n";
    }

    for (unsigned int i = 0; i < reads.size(); i += 2) {
        auto id = reads[i].first;

        int counter = 0;
        for (const auto& path : mPerRead[i]) {
            for (length_t i = 0; i < path.second.size(); i++) {
                // For occurrences shorter than k, use / to show that the nodes
                // do not form a path, but a set of possibilities
                char separationchar =
                    path.second[i].getRange().width() < k_DBG ? '/' : ',';
                f2 << id << "\t" << counter << "\t" << path.first[0];
                for (length_t i = 1; i < path.first.size(); i++) {
                    f2 << separationchar << path.first[i];
                }
                f2 << "\t" << path.second[i].getDistanceFromLeftEnd() << "\t"
                   << path.second[i].getStrain() << "\t"
                   << path.second[i].getRange().getBegin() << "\t"
                   << path.second[i].getRange().width() << "\t"
                   << path.second[i].getDistance() << "\t0\n";
            }
            counter++;
        }

        for (const auto& path : mPerRead[i + 1]) {
            for (length_t i = 0; i < path.second.size(); i++) {
                // For occurrences shorter than k, use / to show that the nodes
                // do not form a path, but a set of possibilities
                char separationchar =
                    path.second[i].getRange().width() < k_DBG ? '/' : ',';
                f2 << id << "\t" << counter << "\t" << path.first[0];
                for (length_t i = 1; i < path.first.size(); i++) {
                    f2 << separationchar << path.first[i];
                }
                f2 << "\t" << path.second[i].getDistanceFromLeftEnd() << "\t"
                   << path.second[i].getStrain() << "\t"
                   << path.second[i].getRange().getBegin() << "\t"
                   << path.second[i].getRange().width() << "\t"
                   << path.second[i].getDistance() << "\t1\n";
            }
            counter++;
        }
    }
}

void writeToOutputSFR(const string& file,
                      const std::vector<std::vector<FMOccSFR>>& mPerRead,
                      const vector<pair<string, string>>& reads,
                      bool& firstChunk, ofstream& f2) {

    if (firstChunk) {
        // Write header
        firstChunk = false;

        // cout << "Writing to output file " << file << " ..." << endl;

        f2 << "Identifier\tSubgraphID\tPath\tDistanceFromLeftEnd\tLength\tED\t"
              "reverseComplement\n";
    }

    for (unsigned int i = 0; i < reads.size(); i += 2) {
        auto id = reads[i].first;

        int counter = 0;
        for (const auto& path : mPerRead[i]) {
            // For occurrences shorter than k, use / to show that the nodes
            // do not form a path, but a set of possibilities
            char separationchar =
                path.getPosition().getTrueDepth() < k_DBG ? '/' : ',';
            vector<uint32_t> nodepath = path.getPosition().getNodePath();
            f2 << id << "\t" << counter << "\t" << nodepath[0];
            for (length_t i = 1; i < nodepath.size(); i++) {
                f2 << separationchar << nodepath[i];
            }
            f2 << "\t" << path.getPosition().getDistanceFromLeftEnd() << "\t"
               << path.getPosition().getTrueDepth() << "\t"
               << path.getDistance() << "\t0\n";
            counter++;
        }

        for (const auto& path : mPerRead[i + 1]) {
            // For occurrences shorter than k, use / to show that the nodes
            // do not form a path, but a set of possibilities
            char separationchar =
                path.getPosition().getTrueDepth() < k_DBG ? '/' : ',';
            vector<uint32_t> nodepath = path.getPosition().getNodePath();
            f2 << id << "\t" << counter << "\t" << nodepath[0];
            for (length_t i = 1; i < nodepath.size(); i++) {
                f2 << separationchar << nodepath[i];
            }
            f2 << "\t" << path.getPosition().getDistanceFromLeftEnd() << "\t"
               << path.getPosition().getTrueDepth() << "\t"
               << path.getDistance() << "\t1\n";
            counter++;
        }
    }
}

double findMedian(vector<length_t> a, int n) {

    // If size of the arr[] is even
    if (n % 2 == 0) {

        // Applying nth_element
        // on n/2th index
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Applying nth_element
        // on (n-1)/2 th index
        nth_element(a.begin(), a.begin() + (n - 1) / 2, a.end());

        // Find the average of value at
        // index N/2 and (N-1)/2
        return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0;
    }

    // If size of the arr[] is odd
    else {

        // Applying nth_element
        // on n/2
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Value at index (N/2)th
        // is the median
        return (double)a[n / 2];
    }
}

template <class T, class positionClass>
double doBenchSFI(T& index, SearchStrategyDBG<T, positionClass>* strategy,
                  string readsFile, length_t ED, string cpSparse,
                  string outputFile) {

    if (outputFile == "") {
        outputFile = readsFile + "_output.tsv";
    }

    size_t totalNodes = 0;
    size_t totalMatrixElements = 0;
    size_t allReportedMatches = 0;
    size_t totalUniqueMatches = 0;
    size_t mappedReads = 0;
    size_t mappedReadsForward = 0;
    size_t mappedReadsBackward = 0;
    size_t totalDBGNodes = 0;
    size_t totalNodePaths = 0;
    size_t allReportedNodePaths = 0;
    chrono::duration<double> elapsedNodepaths =
        std::chrono::duration<double>::zero();
    chrono::duration<double> elapsedSAtoText =
        std::chrono::duration<double>::zero();
    chrono::duration<double> elapsed = std::chrono::duration<double>::zero();

    cout << "Strain-fixed read mapping with " << strategy->getName()
         << " strategy for max distance " << ED << " with "
         << strategy->getPartitioningStrategy() << " partitioning and using "
         << strategy->getDistanceMetric()
         << " distance and using checkpoint sparseness factor " << cpSparse
         << endl;
    cout.precision(2);

    // Read and write input in chunks
    size_t chunkSize = 10000;
    bool firstChunk = true;
    size_t nrOfReads = 0;
    string currentLine = "";

    ifstream ifile(readsFile.c_str());
    ofstream f2;
    f2.open(outputFile);

    bool readWithN = false;

    std::vector<length_t> numberMatchesPerRead;

    while (ifile) {

        // cout << "Reading in reads from " << readsFile << endl;
        vector<pair<string, string>> reads;
        reads.reserve(chunkSize * 2);
        try {
            nrOfReads += getReads(reads, readsFile, ifile, chunkSize,
                                  currentLine, readWithN);
        } catch (const exception& e) {
            string er = e.what();
            er += " Did you provide a valid reads file?";
            throw runtime_error(er);
        }

        auto start = chrono::high_resolution_clock::now();

        vector<std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>>
            matchesPerRead = {};
        matchesPerRead.reserve(reads.size());

        numberMatchesPerRead.reserve(reads.size() +
                                     numberMatchesPerRead.size());

        for (unsigned int i = 0; i < reads.size(); i += 2) {

            const auto& p = reads[i];

            auto originalPos = p.first;
            string read = p.second;
            string revCompl = reads[i + 1].second;

            if ((((i + 2 * nrOfReads - reads.size()) >> 1) - 1) %
                    (8192 / (1 << ED)) ==
                0) {
                cout << "Progress: " << (i + 2 * nrOfReads - reads.size()) / 2
                     << " reads done.\r";
                cout.flush();
            }

            const auto& matches = strategy->matchApproxSFI(read, ED);
            int nr_of_matches = 0;
            for (auto const& p : matches) {
                nr_of_matches += p.second.size();
            }

            totalNodes += index.getNodes();
            totalDBGNodes += index.getDBGNodes();
            totalMatrixElements += index.getMatrixElements();
            allReportedMatches += index.getTotalReported();
            totalUniqueMatches += nr_of_matches;
            mappedReadsForward += !matches.empty();
            totalNodePaths += matches.size();
            allReportedNodePaths += index.getTotalReportedNodePaths();
            elapsedNodepaths += index.getNodePathDuration();
            elapsedSAtoText += index.getSADuration();

            // do the same for the reverse complement
            const auto& matchesRevCompl =
                strategy->matchApproxSFI(revCompl, ED);
            int nr_of_matchesRevCompl = 0;
            for (auto const& p : matchesRevCompl) {
                nr_of_matchesRevCompl += p.second.size();
            }

            totalNodes += index.getNodes();
            totalDBGNodes += index.getDBGNodes();
            totalMatrixElements += index.getMatrixElements();
            allReportedMatches += index.getTotalReported();
            totalUniqueMatches += nr_of_matchesRevCompl;
            mappedReadsBackward += !matchesRevCompl.empty();
            totalNodePaths += matchesRevCompl.size();
            allReportedNodePaths += index.getTotalReportedNodePaths();
            elapsedNodepaths += index.getNodePathDuration();
            elapsedSAtoText += index.getSADuration();

            mappedReads += !(matchesRevCompl.empty() && matches.empty());

            matchesPerRead.push_back(matches);
            matchesPerRead.push_back(matchesRevCompl);
            numberMatchesPerRead.push_back(nr_of_matches +
                                           nr_of_matchesRevCompl);
        }

        auto finish = chrono::high_resolution_clock::now();
        elapsed += finish - start;
        writeToOutputSFI(outputFile, matchesPerRead, reads, firstChunk, f2);
    }

    if (readWithN) {
        cout << "Caution, reads containing an N were ignored.\n";
    }

    f2.close();

    if (ED == 0) {
        allReportedMatches = totalUniqueMatches;
    }

    elapsedSAtoText -= elapsedNodepaths;
    chrono::duration<double> FMIndexElapsed =
        elapsed - elapsedSAtoText - elapsedNodepaths;
    cout << "Progress: " << nrOfReads << "/" << nrOfReads << "\n";
    cout << "Results for " << strategy->getName() << endl;

    cout << "Time for finding the occurrences in the bidirectional FM-index: "
         << fixed << FMIndexElapsed.count() << "s\n";
    cout << "Time for finding the node path: " << fixed
         << elapsedNodepaths.count() << "s\n";
    cout << "Time for finding the occurrences in the reference text along with "
            "the corresponding strain: "
         << fixed << elapsedSAtoText.count() << "s\n";
    cout << "Total duration: " << fixed << elapsed.count() << "s\n";

    cout << "Average no. index nodes: " << totalNodes / (double)(nrOfReads)
         << endl;
    cout << "Total no. index nodes: " << totalNodes << "\n";
    cout << "Average no. unique matches: "
         << totalUniqueMatches / (double)(nrOfReads) << endl;
    cout << "Total no. unique matches: " << totalUniqueMatches << "\n";
    cout << "Average no. reported matches "
         << allReportedMatches / (double)(nrOfReads) << endl;
    cout << "Total no. reported matches: " << allReportedMatches << "\n";
    cout << "Average no. unique node paths "
         << totalNodePaths / (double)(nrOfReads) << endl;
    cout << "Total no. unique node paths: " << totalNodePaths << "\n";
    cout << "Average no. reported node paths: "
         << allReportedNodePaths / (double)(nrOfReads) << endl;
    cout << "Total no. reported node paths: " << allReportedNodePaths << "\n";
    cout << "Mapped reads :" << mappedReads << endl;
    cout << "Median number of unique matches per read "
         << findMedian(numberMatchesPerRead, numberMatchesPerRead.size())
         << endl;
    cout << "Average no. graph nodes: " << totalDBGNodes / (double)(nrOfReads)
         << endl;
    cout << "Total no. graph nodes: " << totalDBGNodes << "\n";

    return elapsed.count();
}

double doBenchSFR(FMIndexDBG<FMPosSFR>& index,
                  SearchStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>* strategy,
                  string readsFile, length_t ED, string cpSparse,
                  string outputFile) {

    StrainFreeMapper mapper(strategy);

    if (outputFile == "") {
        outputFile = readsFile + "_output.tsv";
    }

    size_t totalNodes = 0;
    // size_t totalMatrixElements = 0;
    size_t allReportedMatches = 0;
    size_t totalUniqueMatches = 0;
    size_t mappedReads = 0;
    // size_t mappedReadsForward = 0;
    // size_t mappedReadsBackward = 0;
    size_t totalDBGNodes = 0;
    // size_t totalFilterSpecialCases = 0;
    chrono::duration<double> elapsed = std::chrono::duration<double>::zero();

    cout << "Strain-free read mapping with " << index.getFilteringOption()
         << " filtering, with " << strategy->getName()
         << " strategy for max distance " << ED << " with "
         << strategy->getPartitioningStrategy() << " partitioning, using "
         << strategy->getDistanceMetric()
         << " distance and using checkpoint sparseness factor " << cpSparse
         << endl;
    cout.precision(2);

    size_t chunkSize = 10000;
    bool firstChunk = true;
    size_t nrOfReads = 0;
    string currentLine = "";

    ifstream ifile(readsFile.c_str());
    ofstream f2;
    f2.open(outputFile);

    bool readWithN = false;

    std::vector<length_t> numberMatchesPerRead;

    while (ifile) {

        // cout << "Reading in reads from " << readsFile << endl;
        vector<pair<string, string>> reads;
        reads.reserve(chunkSize * 2);
        try {
            nrOfReads += getReads(reads, readsFile, ifile, chunkSize,
                                  currentLine, readWithN);
        } catch (const exception& e) {
            string er = e.what();
            er += " Did you provide a valid reads file?";
            throw runtime_error(er);
        }
        auto start = chrono::high_resolution_clock::now();

        std::vector<std::vector<FMOccSFR>> matchesPerRead = {};
        matchesPerRead.reserve(reads.size());

        numberMatchesPerRead.reserve(reads.size() +
                                     numberMatchesPerRead.size());

        for (unsigned int i = 0; i < reads.size(); i += 2) {

            const auto& p = reads[i];

            auto originalPos = p.first;
            string read = p.second;
            string revCompl = reads[i + 1].second;

            if ((((i + 2 * nrOfReads - reads.size()) >> 1) - 1) %
                    (8192 / (1 << ED)) ==
                0) {
                cout << "Progress: " << (i + 2 * nrOfReads - reads.size()) / 2
                     << " reads done.\r";
                cout.flush();
            }

            const auto& matches = mapper.matchApproxSFR(read, ED);
            int nr_of_matches = matches.size();

            totalNodes += index.getNodes();
            totalDBGNodes += index.getDBGNodes();
            // totalFilterSpecialCases += index.getFilterSpecialCases();
            // totalMatrixElements += index.getMatrixElements();
            allReportedMatches += index.getTotalReported();
            totalUniqueMatches += nr_of_matches;
            // mappedReadsForward += !matches.empty();

            // do the same for the reverse complement
            const auto& matchesRevCompl = mapper.matchApproxSFR(revCompl, ED);
            int nr_of_matchesRevCompl = matchesRevCompl.size();

            totalNodes += index.getNodes();
            totalDBGNodes += index.getDBGNodes();
            // totalFilterSpecialCases += index.getFilterSpecialCases();
            // totalMatrixElements += index.getMatrixElements();
            allReportedMatches += index.getTotalReported();
            totalUniqueMatches += nr_of_matchesRevCompl;
            // mappedReadsBackward += !matchesRevCompl.empty();

            mappedReads += !(matchesRevCompl.empty() && matches.empty());

            matchesPerRead.push_back(matches);
            matchesPerRead.push_back(matchesRevCompl);
            numberMatchesPerRead.push_back(nr_of_matches +
                                           nr_of_matchesRevCompl);
        }

        auto finish = chrono::high_resolution_clock::now();
        elapsed += finish - start;

        writeToOutputSFR(outputFile, matchesPerRead, reads, firstChunk, f2);
    }

    if (readWithN) {
        cout << "Caution, reads containing an N were ignored.\n";
    }

    f2.close();

    if (ED == 0) {
        allReportedMatches = totalUniqueMatches;
    }

    cout << "Progress: " << nrOfReads << "/" << nrOfReads << "\n";
    cout << "Results for " << strategy->getName() << endl;

    cout << "Total duration: " << fixed << elapsed.count() << "s\n";

    cout << "Average no. index nodes: " << totalNodes / (double)(nrOfReads)
         << endl;
    cout << "Total no. index nodes: " << totalNodes << "\n";
    cout << "Average no. unique node paths: "
         << totalUniqueMatches / (double)(nrOfReads) << endl;
    cout << "Total no. unique node paths: " << totalUniqueMatches << "\n";
    cout << "Average no. reported node paths: "
         << allReportedMatches / (double)(nrOfReads) << endl;
    cout << "Total no. reported node paths: " << allReportedMatches << "\n";
    cout << "Mapped reads :" << mappedReads << endl;
    cout << "Median number of unique node paths per read "
         << findMedian(numberMatchesPerRead, numberMatchesPerRead.size())
         << endl;
    cout << "Average no. graph nodes: " << totalDBGNodes / (double)(nrOfReads)
         << endl;
    cout << "Total no. graph nodes: " << totalDBGNodes << "\n";
    // cout << "Total no. special filter cases: " << totalFilterSpecialCases
    //      << "\n";

    return elapsed.count();
}

template double doBenchSFI<FMIndexDBG<FMPos>, FMPos>(
    FMIndexDBG<FMPos>& index,
    SearchStrategyDBG<FMIndexDBG<FMPos>, FMPos>* strategy, string readsFile,
    length_t ED, string cpSparse, string outputFile);