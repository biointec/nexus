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
#include "searchstrategy.h"
#include <algorithm> // minmax
#include <chrono>    // timings
#include <numeric>   // for summing over vector
#include <sstream>   // used for splitting strings

using namespace std;
vector<string> schemes = {"kuch1",  "kuch2", "kianfar", "manbest",
                          "pigeon", "01*0",  "custom",  "naive"};

int editDistDP(string P, string O, int maxED) {

    int n = (int)P.length();
    int m = (int)O.length();
    string* horizontal = &P;
    string* vertical = &O;
    if (n > m) {
        horizontal = &O;
        vertical = &P;
        int temp = n;
        n = m;
        m = temp;
    }

    // check the dimensions of s1 and s2
    if ((max(m, n) - min(m, n)) > maxED)
        return numeric_limits<int>::max();

    BandMatrix mat(m + 1 + maxED, maxED);

    // fill in the rest of the matrix
    for (int i = 1; i <= m; i++) {
        for (int j = mat.getFirstColumn(i); j <= mat.getLastColumn(i) && j <= n;
             j++) {
            mat.updateMatrix(vertical->at(i - 1) != horizontal->at(j - 1), i,
                             j);
        }
    }
    return mat(m, n);
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
void writeToOutput(const string& file,
                   const vector<vector<TextOccurrence>>& mPerRead,
                   const vector<pair<string, string>>& reads) {

    cout << "Writing to output file " << file << " ..." << endl;
    ofstream f2;
    f2.open(file);

    f2 << "identifier\tposition\tlength\tED\treverseComplement\n";
    for (unsigned int i = 0; i < reads.size(); i += 2) {
        auto id = reads[i].first;

        for (auto m : mPerRead[i]) {
            f2 << id << "\t" << m.getRange().getBegin() << "\t"
               << m.getRange().width() << "\t" << m.getDistance() << "\t0\n";
        }

        for (auto m : mPerRead[i + 1]) {
            f2 << id << "\t" << m.getRange().getBegin() << "\t"
               << m.getRange().width() << "\t" << m.getDistance() << "\t1\n";
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

template <class T>
void doBench(vector<pair<string, string>>& reads, FMIndex<T>& mapper,
             SearchStrategy<FMIndex<FMPos>, FMPos>* strategy, string readsFile,
             length_t ED) {

    size_t totalNodes = 0;
    size_t totalMatrixElements = 0;
    size_t allReportedMatches = 0;
    size_t totalUniqueMatches = 0;
    size_t mappedReads = 0;
    // size_t mappedReadsForward = 0;
    // size_t mappedReadsBackward = 0;

    cout << "Benchmarking with " << strategy->getName()
         << " strategy for max distance " << ED << " with "
         << strategy->getPartitioningStrategy() << " partitioning and using "
         << strategy->getDistanceMetric() << " distance " << endl;
    cout.precision(2);

    vector<vector<TextOccurrence>> matchesPerRead = {};
    matchesPerRead.reserve(reads.size());

    vector<length_t> numberMatchesPerRead;
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

        auto matches = strategy->matchApprox(read, ED);

        totalNodes += mapper.getNodes();
        totalMatrixElements += mapper.getMatrixElements();
        allReportedMatches += mapper.getTotalReported();
        totalUniqueMatches += matches.size();
        // mappedReadsForward += !matches.empty();

        // do the same for the reverse complement
        vector<TextOccurrence> matchesRevCompl =
            strategy->matchApprox(revCompl, ED);
        totalNodes += mapper.getNodes();
        totalMatrixElements += mapper.getMatrixElements();
        allReportedMatches += +mapper.getTotalReported();
        totalUniqueMatches += matchesRevCompl.size();
        // mappedReadsBackward += !matchesRevCompl.empty();

        mappedReads += !(matchesRevCompl.empty() && matches.empty());

        matchesPerRead.push_back(matches);
        matchesPerRead.push_back(matchesRevCompl);
        numberMatchesPerRead.push_back(matches.size() + matchesRevCompl.size());
    }

    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Progress: " << reads.size() << "/" << reads.size() << "\n";
    cout << "Results for " << strategy->getName() << endl;

    cout << "Total duration: " << fixed << elapsed.count() << "s\n";

    cout << "Average no. index nodes: " << totalNodes / (double)(reads.size())
         << endl;
    cout << "Total no. index nodes: " << totalNodes << "\n";
    cout << "Average no. matrix elements written: "
         << totalMatrixElements / (double)(reads.size()) << endl;
    cout << "Total no. matrix elements: " << totalMatrixElements << "\n";
    cout << "Average no. unique matches: "
         << totalUniqueMatches / (double)(reads.size()) << endl;
    cout << "Total no. unique matches: " << totalUniqueMatches << "\n";
    cout << "Average no. reported matches "
         << allReportedMatches / (double)(reads.size()) << endl;
    cout << "Total no. reported matches: " << allReportedMatches << "\n";
    cout << "Mapped reads :" << mappedReads << endl;
    cout << "Median number of occurrences per read "
         << findMedian(numberMatchesPerRead, numberMatchesPerRead.size())
         << endl;

    writeToOutput(readsFile + "_output.tsv", matchesPerRead, reads);
}

void showUsage() {
    cout
        << "Usage: ./columba [options] <basefilename> <k> <readfile.[ext]>\n\n";
    cout << " [options]\n";
    cout << "  -e  --max-ed\t\tmaximum edit distance [default = 0]\n";
    cout << "  -s  --sa-sparseness\tsuffix array sparseness factor "
            "[default = "
            "1]\n";
    cout << "  -p  --partitioning \tAdd flag to do uniform/static/dynamic "
            "partitioning [default = "
            "dynamic]\n";
    cout << "  -m   --metric\tAdd flag to set distance metric "
            "(editnaive/editopt/hamming) [default = "
            "editopt]\n";
    cout << "  -ss --search-scheme\tChoose the search scheme\n  options:\n\t"
         << "kuch1\tKucherov k + 1\n\t"
         << "kuch2\tKucherov k + 2\n\t"
         << "kianfar\t Optimal Kianfar scheme\n\t"
         << "manbest\t Manual best improvement for Kianfar scheme (only for ed "
            "= 4)\n\t"
         << "pigeon\t Pigeonhole scheme\n\t"
         << "01*0\t01*0 search scheme\n\t"
         << "custom\tcustom search scheme, the next parameter should be a path "
            "to the folder containing this search scheme\n\n";

    cout << "[ext]\n"
         << "\tone of the following: fq, fastq, FASTA, fasta, fa\n";

    cout << "Following input files are required:\n";
    cout << "\t<base filename>.txt: input text T\n";
    cout << "\t<base filename>.cct: character counts table\n";
    cout << "\t<base filename>.sa.[saSF]: suffix array sample every [saSF] "
            "elements\n";
    cout << "\t<base filename>.bwt: BWT of T\n";
    cout << "\t<base filename>.brt: Prefix occurrence table of T\n";
    cout << "\t<base filename>.rev.brt: Prefix occurrence table of the "
            "reverse "
            "of T\n";
    cout <<

        "This program aligns short, single end reads to a pan-genome in the\n"
        "form of a linear concatenation. It reports the corresponding\n"
        "coordinates in the original genomes.\n\n\n"

        "Usage: ./columba [options] <basefilename> <readfile.[ext]>\n\n"

        " Following input parameters are required:\n"
        "  <basefilename>      base filename of the input index\n"
        "  <readfile.[ext]>    the file containing the input reads to be\n"
        "                      aligned (single end).\n\n"

        " [ext]\n"
        "  one of the following: fq, fastq, FASTA, fasta, fa\n\n\n"

        " [options]\n"
        "  -e/--max-ed         maximum edit distance [default = 0]\n\n"

        "  -s/--sa-sparseness  suffix array sparseness factor [default = "
        "16]\n\n"

        "  -p/--partitioning   Add flag to do uniform/static/dynamic\n"
        "                      partitioning of the seeds for search schemes.\n"
        "                      Dynamic partitioning cannot be used with\n"
        "                      strain-free matching. [default = dynamic]\n\n"

        "  -m/--metric         Add flag to set distance metric (editnaive/\n"
        "                      editopt/ hamming) [default = editopt]\n\n"

        "  -ss/--search-scheme Choose the search scheme. Options:\n"
        "                       * kuch1    Kucherov k + 1 [default]\n"
        "                       * kuch2    Kucherov k + 2\n"
        "                       * kianfar  Optimal Kianfar scheme\n"
        "                       * manbest  Manual best improvement for "
        "Kianfar\n"
        "                                  scheme (only for ed = 4)\n"
        "                       * pigeon   Pigeonhole scheme\n"
        "                       * 01*0     01*0 search scheme\n"
        "                       * naive    naive backtracking\n"
        "                       * custom   custom search scheme, the next\n"
        "                                  parameter should be a path to the\n"
        "                                  folder containing this "
        "searchscheme\n\n\n"

        " Following input files are required:\n"
        "  <basefilename>.compressed.txt:           compressed version of the\n"
        "                                           input text T\n\n"
        "  <basefilename>.cct:                      character counts table\n\n"
        "  <basefilename>.sa.<saSF>:                sparse suffix array, with\n"
        "                                           suffix array sparseness\n"
        "                                           factor <saSF> elements\n\n"
        "  <basefilename>.sa.bv.<saSF>:             bitvector indicating "
        "which\n"
        "                                           elements of the suffix\n"
        "                                           array are stored.\n\n"
        "  <basefilename>.bwt:                      BWT of T\n\n"
        "  <basefilename>.rev.bwt:                  BWT of the reverse of T\n\n"
        "  <basefilename>.brt:                      Prefix occurrence table of "
        "T\n\n"
        "  <basefilename>.rev.brt:                  Prefix occurrence table "
        "of\n\n\n";
}

int main(int argc, char* argv[]) {

    int requiredArguments = 2; // baseFile of files and file containing reads

    if (argc < requiredArguments) {
        cerr << "Insufficient number of arguments" << endl;
        showUsage();
        return EXIT_FAILURE;
    }
    if (argc == 2 && strcmp("help", argv[1]) == 0) {
        showUsage();
        return EXIT_SUCCESS;
    }

    cout << "Welcome to Columba!\n";

    string saSparse = "1";
    string maxED = "0";
    string searchscheme = "kuch1";
    string customFile = "";

    PartitionStrategy pStrat = DYNAMIC;
    DistanceMetric metric = EDITOPTIMIZED;

    // process optional arguments
    for (int i = 1; i < argc - requiredArguments; i++) {
        const string& arg = argv[i];

        if (arg == "-p" || arg == "--partitioning") {
            if (i + 1 < argc) {
                string s = argv[++i];
                if (s == "uniform") {
                    pStrat = UNIFORM;
                } else if (s == "dynamic") {
                    pStrat = DYNAMIC;
                } else if (s == "static") {
                    pStrat = STATIC;
                } else {
                    throw runtime_error(
                        s + " is not a partitioning option\nOptions are: "
                            "uniform, static, dynamic");
                }

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-s" || arg == "--sa-sparseness") {
            if (i + 1 < argc) {
                saSparse = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-e" || arg == "--max-ed") {
            if (i + 1 < argc) {
                maxED = argv[++i];
            }
        } else if (arg == "-ss" || arg == "--search-scheme") {
            if (i + 1 < argc) {
                searchscheme = argv[++i];
                if (find(schemes.begin(), schemes.end(), searchscheme) ==
                    schemes.end()) {
                    throw runtime_error(searchscheme +
                                        " is not on option as search scheme");
                }
                if (searchscheme == "custom") {
                    if (i + 1 < argc) {
                        customFile = argv[++i];
                    } else {
                        throw runtime_error(
                            "custom search scheme takes a folder as argument");
                    }
                }

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }

        } else if (arg == "-m" || arg == "-metric") {
            if (i + 1 < argc) {
                string s = argv[++i];
                if (s == "editopt") {
                    metric = EDITOPTIMIZED;
                } else if (s == "editnaive") {
                    metric = EDITNAIVE;
                } else if (s == "hamming") {
                    metric = HAMMING;
                } else {
                    throw runtime_error(s +
                                        " is not a metric option\nOptions are: "
                                        "editopt, editnaive, hamming");
                }

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        }

        else {
            cerr << "Unknown argument: " << arg << " is not an option" << endl;
            return EXIT_FAILURE;
        }
    }

    length_t ed = stoi(maxED);
    if (ed < 0 || ed > 5) {
        cerr << ed << " is not allowed as maxED should be in [0, 4]" << endl;

        return EXIT_FAILURE;
    }
    length_t saSF = stoi(saSparse);
    if (saSF == 0 || saSF > 256 || (saSF & (saSF - 1)) != 0) {
        cerr << saSF
             << " is not allowed as sparse factor, should be in 2^[0, 8]"
             << endl;
    }

    if (ed != 4 && searchscheme == "manbest") {
        throw runtime_error("manbest only supports 4 allowed errors");
    }

    string baseFile = argv[argc - 2];
    string readsFile = argv[argc - 1];

    cout << "Reading in reads from " << readsFile << endl;
    vector<pair<string, string>> reads;
    try {
        reads = getReads(readsFile);
    } catch (const exception& e) {
        string er = e.what();
        er += " Did you provide a valid reads file?";
        throw runtime_error(er);
    }
    cout << "Start creation of BWT approximate matcher" << endl;

    FMIndex<FMPos> bwt = FMIndex<FMPos>(baseFile, saSF);

    SearchStrategy<FMIndex<FMPos>, FMPos>* strategy;
    if (searchscheme == "kuch1") {
        strategy =
            new KucherovKplus1<FMIndex<FMPos>, FMPos>(bwt, pStrat, metric);
    } else if (searchscheme == "kuch2") {
        strategy =
            new KucherovKplus2<FMIndex<FMPos>, FMPos>(bwt, pStrat, metric);
    } else if (searchscheme == "kianfar") {
        strategy =
            new OptimalKianfar<FMIndex<FMPos>, FMPos>(bwt, pStrat, metric);
    } else if (searchscheme == "manbest") {
        strategy =
            new ManBestStrategy<FMIndex<FMPos>, FMPos>(bwt, pStrat, metric);
    } else if (searchscheme == "01*0") {
        strategy = new O1StarSearchStrategy<FMIndex<FMPos>, FMPos>(bwt, pStrat,
                                                                   metric);
    } else if (searchscheme == "pigeon") {
        strategy = new PigeonHoleSearchStrategy<FMIndex<FMPos>, FMPos>(
            bwt, pStrat, metric);
    } else if (searchscheme == "custom") {
        strategy = new CustomSearchStrategy<FMIndex<FMPos>, FMPos>(
            bwt, customFile, pStrat, metric);
    } else if (searchscheme == "naive") {
        strategy = new NaiveBackTrackingStrategy<FMIndex<FMPos>, FMPos>(
            bwt, pStrat, metric);
    } else {
        // should not get here
        throw runtime_error(searchscheme +
                            " is not on option as search scheme");
    }
    doBench(reads, bwt, strategy, readsFile, ed);
    delete strategy;
    cout << "Bye...\n";
}
