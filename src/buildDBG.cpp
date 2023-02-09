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

#include <sstream>

using namespace std;

typedef uint64_t length_t;

void showUsage() {
    cout
        << "This program constructs a new implicit pan-genome de Bruijn\n"
           "graph, along with the underlying bidirectional FM-index.\n\n\n"

           "Usage: ./nexusBuild <base filename> <k_list>\n\n"

           " Following input parameters are required:\n"
           "  <base filename>       base filename of the input text\n"
           "  <k_list>              the de Bruijn parameter: a "
           "comma-separated\n"
           "                        list of integers is required (e.g.,\n"
           "                        20,21,23)\n\n"

           " Following input files are required:\n"
           "  <base filename>.txt:  input text with all genomes readily\n"
           "                        concatenated, containing only the "
           "following\n"
           "                        characters: A, C, G, T, % and $ (at the\n"
           "                        very end). No newlines are allowed.\n\n\n"

           " [options]\n"
           "  --skip                Skip the building process of the data\n"
           "                        structures that are independent of the de\n"
           "                        Bruijn k parameter (i.e., the "
           "bidirectional\n"
           "                        FM-index). These data structures must be\n"
           "                        available in the directory.\n\n"
           "  -s/--sa-sparseness    Suffix array sparseness factors to be\n"
           "                        used. This option can be repeated "
           "multiple\n"
           "                        times for multiple versions of the suffix\n"
           "                        array. This option takes values in {1, 2,\n"
           "                        4, 8, 16, 32, 64, 128, 256}. Use \"all\"\n"
           "                        to use all options. [default = 16]\n\n"
           "  -c/--cp-sparseness    Sparseness factor that indicates how many\n"
           "                        checkpoints must be stored to identify\n"
           "                        nodes. This option can be repeated\n"
           "                        multiple times for multiple versions of\n"
           "                        the checkpoint sparseness. Use \"none\" "
           "to\n"
           "                        use no checkpoints. [default = 128]\n\n"
           "  -p  --progress        Report extra progress updates\n\n\n";
}

bool parseArguments(int argc, char* argv[], string& baseFN, vector<uint>& k,
                    vector<int>& saSF, vector<int>& cpSF, bool& progress,
                    bool& skip) {
    const int reqArguments = 2;
    progress = false;

    if (argc == 2) {
        string firstArg(argv[1]);
        if (firstArg.find("help") != std::string::npos) {
            showUsage();
            return EXIT_SUCCESS;
        }
    }

    if (argc <= reqArguments) {
        cerr << "Fatal error: insufficient number of arguments.\n" << endl;
        return false;
    }
    // process optional arguments
    for (int i = 1; i < argc - reqArguments; i++) {
        string arg(argv[i]);

        // process options
        if (((arg == "-s") || (arg == "--sa-sparseness"))) {
            string arg2(argv[i + 1]);
            if (arg2 == "all") {
                saSF.insert(saSF.end(), {256, 128, 64, 32, 16, 8, 4, 2, 1});
            } else {
                int SF = atoi(argv[i + 1]);

                if ((SF != 1) && (SF != 2) && (SF != 4) && (SF != 8) &&
                    (SF != 16) && (SF != 32) && (SF != 64) && (SF != 128) &&
                    (SF != 256)) {
                    cerr << "Suffix array sparseness should take values in {1, "
                            "2, 4, 8, "
                            "16, 32, 64, 128, 256}"
                         << endl;
                    return false;
                }

                saSF.push_back(SF);
            }

            i++;
        } else if (((arg == "-c") || (arg == "--cp-sparseness"))) {
            string arg2(argv[i + 1]);
            if (arg2 == "none") {
                cpSF.push_back(INT32_MAX);
            } else {
                int SF = atoi(argv[i + 1]);
                double logSF = log2(SF);
                double value;

                if (modf(logSF, &value) != 0.0) {
                    cerr << "Checkpoint sparseness should be a power of 2."
                         << endl;
                    return false;
                }

                cpSF.push_back(SF);
            }
            i++;
        } else if (((arg == "-p") || (arg == "--progress"))) {
            progress = true;
        } else if (((arg == "--skip"))) {
            skip = true;
        } else {
            cerr << "Unknown argument: " << argv[i] << endl;
            return false;
        }
    }

    if (saSF.empty()) {
        saSF.push_back(16);
    }

    if (cpSF.empty()) {
        cpSF.push_back(128);
    }

    std::sort(saSF.begin(), saSF.end(), std::greater<int>());
    saSF.erase(unique(saSF.begin(), saSF.end()), saSF.end());

    std::sort(cpSF.begin(), cpSF.end(), std::greater<int>());
    cpSF.erase(unique(cpSF.begin(), cpSF.end()), cpSF.end());

    baseFN = argv[argc - 2];
    string k_list = argv[argc - 1];

    std::stringstream ss(k_list);
    string tmp;
    while (getline(ss, tmp, ',')) {
        k.push_back(stoi(tmp));
    }
    return true;
}

int main(int argc, char* argv[]) {
    string baseFilename;
    vector<uint> k;
    vector<int> saSF = {};
    vector<int> cpSF = {};
    bool progress;
    bool skip = false;

    if (!parseArguments(argc, argv, baseFilename, k, saSF, cpSF, progress,
                        skip)) {
        showUsage();
        return EXIT_FAILURE;
    }

    cout << "Welcome to Nexus!\n";
    cout << "Alphabet size is " << ALPHABET - 1 << " + 1\n";
    // cout << "k is " << k << "\n";

    try {
        // cout << "Start creation of BWT approximate matcher" << endl;

        FMIndexDBG<FMPos>::buildFMIndexDBG(baseFilename, k, saSF, cpSF,
                                           progress, skip);
    } catch (const std::exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    cout << "Exiting... bye!" << endl;
    return EXIT_SUCCESS;
}