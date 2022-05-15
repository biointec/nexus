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
#include "searchstrategy.h"
#include "strainfreemapper.h"

#include <sstream>

/**
 * @brief Show the usage in terminal
 *
 */
void showUsage() {
    cout << "Usage: ./visualizePath [options] basefilename path\n\n";
    cout << " [path] should be a comma-separated list of node identifiers\n\n";
    cout << " [Visualization options]\n";
    cout << " -d   --visualization-depth\t\tDepth of the visualized "
            "neighborhood around the paths of interest [default = 3]\n";
    cout << " -o   --output-files\t\tPrefix of the output files that will be "
            "created during the visualization process [default = "
            "basefilename]\n";
    cout << " -s   --sa-sparseness\tsuffix array sparseness factor "
            "[default = "
            "1]\n";
    cout << " -c   --cp-sparseness\tsparseness factor that indicates "
            "how many checkpoints must be stored to identify nodes. Use "
            "\"none\" to use no checkpoints. Choose a value that was also used "
            "during the building process."
            "[default = 128]\n\n";

    cout << "Following input files are required:\n";
    cout << "\t<base filename>.txt: input text T\n";
    cout << "\t<base filename>.cct: character counts table\n";
    cout << "\t<base filename>.sa.[saSF]: sparse suffix array, with suffix "
            "array sparseness factor [saSF] "
            "elements\n";
    cout << "\t<base filename>.sa.bv.[saSF]: bitvector indicating which "
            "elements of the suffix array are stored.\n";
    cout << "\t<base filename>.bwt: BWT of T\n";
    cout << "\t<base filename>.rev.bwt: BWT of the reverse of T\n";
    cout << "\t<base filename>.brt: Prefix occurrence table of T\n";
    cout << "\t<base filename>.rev.brt: Prefix occurrence table of the "
            "reverse "
            "of T\n";
    cout << "\t<base filename>.DBG: variable k and the compressed de "
            "Bruijn graph.\n";
    cout << "\t<base filename>.B.left: bitvector B_left for the compressed de "
            "Bruijn graph.\n";
    cout << "\t<base filename>.B.right.[cpSF]: bitvector B_right for the "
            "compressed "
            "de "
            "Bruijn graph, with checkpoint sparseness factor [cpSF].\n";
    cout << "\t<base filename>.B.right.full.[cpSF]: bitvector B_right_full for "
            "the "
            "compressed de Bruijn graph, with checkpoint sparseness factor "
            "[cpSF].\n";
    cout << "\t<base filename>.left.map: node identifier mapping corresponding "
            "to B_left.\n";
    cout << "\t<base filename>.right.map.[cpSF]: node identifier mapping "
            "corresponding "
            "to B_right, with checkpoint sparseness factor [cpSF].\n";
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

    cout << "Welcome to Nexus!\n";

    string visDepthString = "3";
    string outputFile = "";
    string cpSparse = "128";
    string saSparse = "1";

    // process optional arguments
    for (int i = 1; i < argc - requiredArguments; i++) {
        const string& arg = argv[i];

        if (arg == "-d" || arg == "--visualization-depth") {
            if (i + 1 < argc) {
                visDepthString = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-o" || arg == "--output-files") {
            if (i + 1 < argc) {
                outputFile = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-c" || arg == "--cp-sparseness") {
            if (i + 1 < argc) {
                cpSparse = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-s" || arg == "--sa-sparseness") {
            if (i + 1 < argc) {
                saSparse = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else {
            cerr << "Unknown argument: " << arg << " is not an option" << endl;
            return false;
        }
    }
    length_t visDepth = stoi(visDepthString);
    length_t cpSF;
    if (cpSparse == "none") {
        cpSF = INT32_MAX;
    } else {
        cpSF = stoi(cpSparse);
        double logcpSF = log2(cpSF);
        double value;
        if (modf(logcpSF, &value) != 0.0) {
            cerr << "Checkpoint sparseness should be a power of 2." << endl;
            return EXIT_FAILURE;
        }
    }
    length_t saSF = stoi(saSparse);
    if (saSF == 0 || saSF > 256 || (saSF & (saSF - 1)) != 0) {
        cerr << saSF
             << " is not allowed as sparse factor, should be in 2^[0, 8]"
             << endl;
        return EXIT_FAILURE;
    }

    string baseFile = argv[argc - 2];
    string pathString = argv[argc - 1];

    if (outputFile == "") {
        outputFile = baseFile;
    }

    cout << "Start creation of BWT approximate matcher on graphs" << endl;

    FMIndexDBG<FMPosSFR> bwt(baseFile, saSF, cpSF);

    std::vector<uint32_t> path;

    try {
        std::stringstream ss(pathString);

        for (size_t i; ss >> i;) {
            path.push_back(i);
            if (ss.peek() == ',')
                ss.ignore();
        }
    } catch (const std::exception& e) {
        std::cerr << "Something went wrong whilst parsing the path." << '\n';
    }

    bwt.getText();

    bwt.visualizeSubgraph(path, visDepth, outputFile);

    cout << "Bye...\n";
}
