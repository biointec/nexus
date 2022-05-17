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

template <class T>
void printMatches(vector<TextOccurrence> matches, string text, bool printLine,
                  string duration, FMIndex<T>& index, string name) {

    cout << endl;

    cout << name << ":\tduration: " << duration
         << "µs\t nodes visited: " << index.getNodes()
         << "\t matrix elements written: " << index.getMatrixElements()
         << "\t startpositions reported: " << index.getTotalReported()
         << " #matches: " << matches.size() << endl;

    for (auto match : matches) {
        cout << "Found match at position " << match.getRange().getBegin()
             << " with ED " << match.getDistance() << endl;

        cout << "\tCorresponding substring:\t"
             << text.substr(match.getRange().getBegin(),
                            match.getRange().getEnd() -
                                match.getRange().getBegin())
             << endl;
    }
}

string getFileExt(const string& s) {

    size_t i = s.rfind('.', s.length());
    if (i != string::npos) {
        return (s.substr(i + 1, s.length() - i));
    }

    return ("");
}

vector<pair<string, string>> getReads(const string& file) {
    vector<pair<string, string>> reads;
    reads.reserve(200000);

    const auto& extension = getFileExt(file);

    bool fasta =
        (extension == "FASTA") || (extension == "fasta") || (extension == "fa");
    bool fastq = (extension == "fq") || (extension == "fastq");

    ifstream ifile(file.c_str());
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
        getline(ifile, line);

        while (getline(ifile, line)) {
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
        // fasta file
        string read = "";
        string p = "";
        string line;
        while (getline(ifile, line)) {
            if (!line.empty() && line[0] == '>') {

                if (!read.empty()) {

                    reads.push_back(make_pair(p, read));
                    reads.push_back(
                        make_pair(p, Nucleotide::getRevCompl(read)));
                    read.clear();
                }

                p = (line.substr(1));

            } else {
                read += line;
            }
        }
        if (!read.empty()) {

            reads.push_back(make_pair(p, read));
            reads.push_back(make_pair(p, Nucleotide::getRevCompl(read)));
            read.clear();
        }
    } else {
        // fastQ
        string read = "";
        string id = "";
        string line;
        bool readLine = false;
        while (getline(ifile, line)) {
            if (!line.empty() && line[0] == '@') {
                if (!read.empty()) {

                    reads.push_back(make_pair(id, read));
                    reads.push_back(
                        make_pair(id, Nucleotide::getRevCompl(read)));
                    read.clear();
                }
                id = (line.substr(1));
                readLine = true;
            } else if (readLine) {
                read = line;
                readLine = false;
            }
        }
        if (!read.empty()) {

            reads.push_back(make_pair(id, read));
            reads.push_back(make_pair(id, Nucleotide::getRevCompl(read)));
            read.clear();
        }
    }

    return reads;
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
    const vector<pair<string, string>>& reads) {

    cout << "Writing to output file " << file << " ..." << endl;
    ofstream f2;
    f2.open(file);

    f2 << "Identifier\tSubgraphID\tPath\tStrain\tPosition\tLength\tED\treverseC"
          "omplement\n";
    for (unsigned int i = 0; i < reads.size(); i += 2) {
        auto id = reads[i].first;

        int counter = 0;
        for (const auto& path : mPerRead[i]) {
            for (length_t i = 0; i < path.second.size(); i++) {
                f2 << id << "\t" << counter << "\t" << path.first[0];
                for (length_t i = 1; i < path.first.size(); i++) {
                    f2 << "," << path.first[i];
                }
                f2 << "\t" << path.second[i].getStrain() << "\t"
                   << path.second[i].getRange().getBegin() << "\t"
                   << path.second[i].getRange().width() << "\t"
                   << path.second[i].getDistance() << "\t0\n";
            }
            counter++;
        }

        for (const auto& path : mPerRead[i + 1]) {
            for (length_t i = 0; i < path.second.size(); i++) {
                f2 << id << "\t" << counter << "\t" << path.first[0];
                for (length_t i = 1; i < path.first.size(); i++) {
                    f2 << "," << path.first[i];
                }
                f2 << "\t" << path.second[i].getStrain() << "\t"
                   << path.second[i].getRange().getBegin() << "\t"
                   << path.second[i].getRange().width() << "\t"
                   << path.second[i].getDistance() << "\t1\n";
            }
            counter++;
        }
    }
    f2.close();
}

void writeToOutputSFR(const string& file,
                      const std::vector<std::vector<FMOccSFR>>& mPerRead,
                      const vector<pair<string, string>>& reads,
                      FMIndexDBG<FMPosSFR>& index) {

    cout << "Writing to output file " << file << " ..." << endl;
    ofstream f2;
    f2.open(file);

    f2 << "Identifier\tSubgraphID\tPath\tDistanceFromLeftEnd\tLength\tED\trever"
          "seC"
          "omplement\n";
    for (unsigned int i = 0; i < reads.size(); i += 2) {
        auto id = reads[i].first;

        int counter = 0;
        for (const auto& path : mPerRead[i]) {
            vector<uint32_t> nodepath = path.getPosition().getNodePath();
            f2 << id << "\t" << counter << "\t" << nodepath[0];
            for (length_t i = 1; i < nodepath.size(); i++) {
                f2 << "," << nodepath[i];
            }
            f2 << "\t" << path.getPosition().getDistanceFromLeftEnd() << "\t"
               << path.getPosition().getTrueDepth() << "\t"
               << path.getDistance() << "\t0\n";
            counter++;
        }

        for (const auto& path : mPerRead[i + 1]) {
            vector<uint32_t> nodepath = path.getPosition().getNodePath();
            f2 << id << "\t" << counter << "\t" << nodepath[0];
            for (length_t i = 1; i < nodepath.size(); i++) {
                f2 << "," << nodepath[i];
            }
            // int id, l;
            // index.findID(path.getPosition().getRanges().getRangeSA().getBegin(),
            //              id, l);
            f2 << "\t" << path.getPosition().getDistanceFromLeftEnd() << "\t"
               << path.getPosition().getTrueDepth() << "\t"
               << path.getDistance() << "\t1\n";
            counter++;
        }
    }
    f2.close();
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
double doBenchSFI(vector<pair<string, string>>& reads, T& index,
                  SearchStrategyDBG<T, positionClass>* strategy,
                  string readsFile, length_t ED, string cpSparse,
                  string outputFile) {

    if (outputFile == "") {
        outputFile = readsFile + "_output.txt";
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

    cout << "Strain-fixed read mapping with " << strategy->getName()
         << " strategy for max distance " << ED << " with "
         << strategy->getPartitioningStrategy() << " partitioning and using "
         << strategy->getDistanceMetric()
         << " distance and using checkpoint sparseness factor " << cpSparse
         << endl;
    cout.precision(2);

    vector<std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>>
        matchesPerRead = {};
    matchesPerRead.reserve(reads.size());

    std::vector<length_t> numberMatchesPerRead;
    numberMatchesPerRead.reserve(reads.size());

    auto start = chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < reads.size(); i += 2) {

        const auto& p = reads[i];

        auto originalPos = p.first;
        string read = p.second;
        string revCompl = reads[i + 1].second;

        if (((i >> 1) - 1) % (8192 / (1 << ED)) == 0) {
            cout << "Progress: " << i / 2 << "/" << reads.size() / 2 << "\r";
            cout.flush();
        }

        std::map<std::vector<uint32_t>, std::vector<TextOccurrenceSFI>>
            matches = strategy->matchApproxSFI(read, ED);
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

        // do the same for the reverse complement
        auto matchesRevCompl = strategy->matchApproxSFI(revCompl, ED);
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

        mappedReads += !(matchesRevCompl.empty() && matches.empty());

        matchesPerRead.push_back(matches);
        matchesPerRead.push_back(matchesRevCompl);
        numberMatchesPerRead.push_back(nr_of_matches + nr_of_matchesRevCompl);
    }

    if (ED == 0) {
        allReportedMatches = totalUniqueMatches;
    }

    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Progress: " << reads.size() << "/" << reads.size() << "\n";
    cout << "Results for " << strategy->getName() << endl;

    cout << "Total duration: " << fixed << elapsed.count() << "s\n";

    cout << "Average no. index nodes: " << totalNodes / (double)(reads.size())
         << endl;
    cout << "Total no. index nodes: " << totalNodes << "\n";
    cout << "Average no. unique matches: "
         << totalUniqueMatches / (double)(reads.size()) << endl;
    cout << "Total no. unique matches: " << totalUniqueMatches << "\n";
    cout << "Average no. reported matches "
         << allReportedMatches / (double)(reads.size()) << endl;
    cout << "Total no. reported matches: " << allReportedMatches << "\n";
    cout << "Average no. unique node paths "
         << totalNodePaths / (double)(reads.size()) << endl;
    cout << "Total no. unique node paths: " << totalNodePaths << "\n";
    cout << "Average no. reported node paths: "
         << allReportedNodePaths / (double)(reads.size()) << endl;
    cout << "Total no. reported node paths: " << allReportedNodePaths << "\n";
    cout << "Mapped reads :" << mappedReads << endl;
    cout << "Median number of unique matches per read "
         << findMedian(numberMatchesPerRead, numberMatchesPerRead.size())
         << endl;
    cout << "Average no. graph nodes: "
         << totalDBGNodes / (double)(reads.size()) << endl;
    cout << "Total no. graph nodes: " << totalDBGNodes << "\n";

    writeToOutputSFI(outputFile, matchesPerRead, reads);
    return elapsed.count();
}

double doBenchSFR(vector<pair<string, string>>& reads,
                  FMIndexDBG<FMPosSFR>& index,
                  SearchStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>* strategy,
                  string readsFile, length_t ED, string cpSparse,
                  string outputFile) {

    StrainFreeMapper mapper(strategy);

    if (outputFile == "") {
        outputFile = readsFile + "_output.txt";
    }

    size_t totalNodes = 0;
    size_t totalMatrixElements = 0;
    size_t allReportedMatches = 0;
    size_t totalUniqueMatches = 0;
    size_t mappedReads = 0;
    size_t mappedReadsForward = 0;
    size_t mappedReadsBackward = 0;
    size_t totalDBGNodes = 0;
    size_t totalFilterSpecialCases = 0;

    cout << "Strain-free read mapping with " << index.getFilteringOption()
         << " filtering, with " << strategy->getName()
         << " strategy for max distance " << ED << " with "
         << strategy->getPartitioningStrategy() << " partitioning, using "
         << strategy->getDistanceMetric()
         << " distance and using checkpoint sparseness factor " << cpSparse
         << endl;
    cout.precision(2);

    std::vector<std::vector<FMOccSFR>> matchesPerRead = {};
    matchesPerRead.reserve(reads.size());

    std::vector<length_t> numberMatchesPerRead;
    numberMatchesPerRead.reserve(reads.size());

    auto start = chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < reads.size(); i += 2) {

        const auto& p = reads[i];

        auto originalPos = p.first;
        string read = p.second;
        string revCompl = reads[i + 1].second;

        if (((i >> 1) - 1) % (8192 / (1 << ED)) == 0) {
            cout << "Progress: " << i / 2 << "/" << reads.size() / 2 << "\r";
            cout.flush();
        }

        std::vector<FMOccSFR> matches = mapper.matchApproxSFR(read, ED);
        int nr_of_matches = matches.size();

        totalNodes += index.getNodes();
        totalDBGNodes += index.getDBGNodes();
        totalFilterSpecialCases += index.getFilterSpecialCases();
        totalMatrixElements += index.getMatrixElements();
        allReportedMatches += index.getTotalReported();
        totalUniqueMatches += nr_of_matches;
        mappedReadsForward += !matches.empty();

        // do the same for the reverse complement
        auto matchesRevCompl = mapper.matchApproxSFR(revCompl, ED);
        int nr_of_matchesRevCompl = matchesRevCompl.size();

        totalNodes += index.getNodes();
        totalDBGNodes += index.getDBGNodes();
        totalFilterSpecialCases += index.getFilterSpecialCases();
        totalMatrixElements += index.getMatrixElements();
        allReportedMatches += index.getTotalReported();
        totalUniqueMatches += nr_of_matchesRevCompl;
        mappedReadsBackward += !matchesRevCompl.empty();

        mappedReads += !(matchesRevCompl.empty() && matches.empty());

        matchesPerRead.push_back(matches);
        matchesPerRead.push_back(matchesRevCompl);
        numberMatchesPerRead.push_back(nr_of_matches + nr_of_matchesRevCompl);
    }

    if (ED == 0) {
        allReportedMatches = totalUniqueMatches;
    }

    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Progress: " << reads.size() << "/" << reads.size() << "\n";
    cout << "Results for " << strategy->getName() << endl;

    cout << "Total duration: " << fixed << elapsed.count() << "s\n";

    cout << "Average no. index nodes: " << totalNodes / (double)(reads.size())
         << endl;
    cout << "Total no. index nodes: " << totalNodes << "\n";
    cout << "Average no. unique node paths: "
         << totalUniqueMatches / (double)(reads.size()) << endl;
    cout << "Total no. unique node paths: " << totalUniqueMatches << "\n";
    cout << "Average no. reported node paths: "
         << allReportedMatches / (double)(reads.size()) << endl;
    cout << "Total no. reported node paths: " << allReportedMatches << "\n";
    cout << "Mapped reads :" << mappedReads << endl;
    cout << "Median number of unique node paths per read "
         << findMedian(numberMatchesPerRead, numberMatchesPerRead.size())
         << endl;
    cout << "Average no. graph nodes: "
         << totalDBGNodes / (double)(reads.size()) << endl;
    cout << "Total no. graph nodes: " << totalDBGNodes << "\n";
    // cout << "Total no. special filter cases: " << totalFilterSpecialCases
    //      << "\n";

    writeToOutputSFR(outputFile, matchesPerRead, reads, index);
    return elapsed.count();
}

template double doBenchSFI<FMIndexDBG<FMPos>, FMPos>(
    vector<pair<string, string>>& reads, FMIndexDBG<FMPos>& index,
    SearchStrategyDBG<FMIndexDBG<FMPos>, FMPos>* strategy, string readsFile,
    length_t ED, string cpSparse, string outputFile);