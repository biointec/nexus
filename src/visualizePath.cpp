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
    cout <<

        "This program visualizes a node, a node path or a set of nodes of\n"
        "interest.\n\n\n"

        "Usage: ./visualizePath [options] <basefilename> <k> <path>\n\n"

        " Following input parameters are required:\n"
        "  <basefilename>      base filename of the input index\n"
        "  <k>                 the de Bruijn parameter of the index\n"
        "  <path>              a comma-separated list of node identifiers\n"
        "                      (e.g., 1,9,20)\n\n\n"

        " [options]\n"
        "  -e/--max-ed         maximum edit distance [default = 0]\n\n"

        "  -s/--sa-sparseness  suffix array sparseness factor [default = "
        "16]\n\n"

        "  -c/--cp-sparseness  sparseness factor that indicates how many\n"
        "                      checkpoints must be stored to identify nodes.\n"
        "                      Use \"none\" to use no checkpoints. Choose a\n"
        "                      value that was also used during the building\n"
        "                      process. [default = 128]\n\n"
        "  -d/--depth          Depth of the visualized neighborhood around "
        "the\n"
        "                      paths of interest [default = 3]\n\n"
        "  -b/--bundle-edges   Bundle edges stemming from different strains\n"
        "                      together. Recommended when many strains are\n"
        "                      present [default = false]\n\n"
        "  -o/--output-files   Prefix of the output files that will be "
        "created\n"
        "                      during the visualization process [default =\n"
        "                      basefilename]\n\n\n"

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
        "of\n"
        "                                           the reverse of T\n\n"
        "  <basefilename>.DBG.k<k>:                 the compressed de Bruijn\n"
        "                                           graph for the requested "
        "de\n"
        "                                           Bruijn parameter\n\n"
        "  <basefilename>.B.right.k<k>.cp<cpSF>:    first bitvector of the\n"
        "                                           implicit representation "
        "for\n"
        "                                           the requested de Bruijn\n"
        "                                           parameter, with "
        "checkpoint\n"
        "                                           sparseness factor "
        "<cpSF>\n\n"
        "  <basefilename>.B.left.k<k>:              second bitvector of the\n"
        "                                           implicit representation\n"
        "                                           for the requested de\n"
        "                                           Bruijn parameter\n\n"
        "  <basefilename>.right.map.k<k>.cp<cpSF>:  node identifier mapping\n"
        "                                           corresponding to the "
        "first\n"
        "                                           bitvector, with "
        "checkpoint\n"
        "                                           sparseness factor "
        "<cpSF>\n\n"
        "  <basefilename>.left.map.k<k>:            node identifier mapping\n"
        "                                           corresponding to the\n"
        "                                           second bitvector\n\n\n";
}

int main(int argc, char* argv[]) {

    int requiredArguments = 3; // baseFile of files, k and the node path

    if (argc == 2) {
        string firstArg(argv[1]);
        if (firstArg.find("help") != std::string::npos) {
            showUsage();
            return EXIT_SUCCESS;
        }
    }

    if (argc < requiredArguments) {
        cerr << "Insufficient number of arguments.\n" << endl;
        showUsage();
        return EXIT_FAILURE;
    }

    cout << "Welcome to Nexus!\n";

    string visDepthString = "3";
    string outputFile = "";
    string cpSparse = "128";
    string saSparse = "16";
    bool separateEdges = true;

    // process optional arguments
    for (int i = 1; i < argc - requiredArguments; i++) {
        const string& arg = argv[i];

        if (arg == "-d" || arg == "--depth") {
            if (i + 1 < argc) {
                visDepthString = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-b" || arg == "--bundle-edges") {
            separateEdges = false;
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
            return EXIT_FAILURE;
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

    string baseFile = argv[argc - 3];
    uint k = atoi(argv[argc - 2]);
    string pathString = argv[argc - 1];

    if (outputFile == "") {
        outputFile = baseFile;
    }

    cout << "Start creation of BWT approximate matcher on graphs" << endl;

    FMIndexDBG<FMPosSFR> bwt(baseFile, saSF, cpSF, k);

    std::vector<uint32_t> path;

    try {
        std::stringstream ss(pathString);

        for (size_t i; ss >> i;) {
            path.push_back(i);
            if (ss.peek() == ',')
                ss.ignore();
        }
    } catch (const std::exception& e) {
        std::cerr << "Something went wrong whilst parsing the node path."
                  << '\n';
    }

    bwt.getText();

    bwt.visualizeSubgraph(path, visDepth, outputFile, separateEdges);

    cout << "Bye...\n";
}
