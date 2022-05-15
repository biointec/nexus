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

using namespace std;

typedef uint64_t length_t;

void showUsage() {
    cout << "Usage: ./nexusBuild <base filename> <k>\n\n";

    cout << "Following input files are required:\n";
    cout << "\t<base filename>.txt: input text T\n\n";

    cout << "Following parameters are required:\n";
    cout << "\t<k> is the de Bruijn parameter\n\n";

    cout << " [options]\n";
    cout << "  -s  --sa-sparseness\tSuffix array sparseness factors to be "
            "used. This option can be repeated multiple times for multiple "
            "versions of the suffix array. This option takes values in {1, 2, "
            "4, 8, 16, 32, 64, 128, 256}. Use \"all\" to use all options. "
            "[default = 1]\n";
    cout << "  -c  --cp-sparseness\tSparseness factor that indicates "
            "how many checkpoints must be stored to identify nodes. This "
            "option can be repeated multiple times for multiple "
            "versions of the checkpoint sparseness. Use \"none\" to use no "
            "checkpoints. "
            "[default = 128]\n";
    cout << "  -p  --progress\tReport extra progress updates\n";
}

bool parseArguments(int argc, char* argv[], string& baseFN, uint& k,
                    vector<int>& saSF, vector<int>& cpSF, bool& progress) {
    const int reqArguments = 2;
    progress = false;
    if (argc == 2 && strcmp("help", argv[1]) == 0) {
        return false;
    }
    if (argc <= reqArguments) {
        cerr << "Fatal error: insufficient number of arguments" << endl;
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
        } else {
            cerr << "Unknown argument: " << argv[i] << endl;
            return false;
        }
    }

    if (saSF.empty()) {
        saSF.push_back(1);
    }

    if (cpSF.empty()) {
        cpSF.push_back(128);
    }

    std::sort(saSF.begin(), saSF.end(), std::greater<int>());
    saSF.erase(unique(saSF.begin(), saSF.end()), saSF.end());

    std::sort(cpSF.begin(), cpSF.end(), std::greater<int>());
    cpSF.erase(unique(cpSF.begin(), cpSF.end()), cpSF.end());

    baseFN = argv[argc - 2];
    k = atoi(argv[argc - 1]);
    return true;
}

int main(int argc, char* argv[]) {
    string baseFilename;
    uint k;
    vector<int> saSF = {};
    vector<int> cpSF = {};
    bool progress;

    if (!parseArguments(argc, argv, baseFilename, k, saSF, cpSF, progress)) {
        showUsage();
        return EXIT_FAILURE;
    }

    cout << "Welcome to Nexus!\n";
    cout << "Alphabet size is " << ALPHABET - 1 << " + 1\n";
    cout << "k is " << k << "\n";

    try {
        // cout << "Start creation of BWT approximate matcher" << endl;

        FMIndexDBG<FMPos>::buildFMIndexDBG(baseFilename, k, saSF, cpSF,
                                           progress);
    } catch (const std::exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    cout << "Exiting... bye!" << endl;
    return EXIT_SUCCESS;
}