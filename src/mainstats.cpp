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

// #include "benchmarking.h"
// #include "searchstrategy.h"
#include "strainfreemapper.h"

using namespace std;

/**
 * @brief Show the usage in terminal
 *
 */
void showUsage() {
    cout <<

        "This program reports some statistics regarding the pan-genome graph\n"
        "topology.\n\n\n"

        "Usage: ./nexusStats [options] <basefilename> <k> \n\n"

        " Following input parameters are required:\n"
        "  <basefilename>      base filename of the input index\n"
        "  <k>                 the de Bruijn parameter of the index\n\n\n"

        " [options]\n"

        "  -s/--sa-sparseness  suffix array sparseness factor [default =\n"
        "                      256 to limit memory usage]\n\n"

        "  -c/--cp-sparseness  sparseness factor that indicates how many\n"
        "                      checkpoints must be stored to identify nodes.\n"
        "                      Use \"none\" to use no checkpoints. Choose a\n"
        "                      value that was also used during the building\n"
        "                      process. [default = 128]\n\n\n"

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

    int requiredArguments = 2; // baseFile of files and k

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

    string saSparse = "256";
    string cpSparse = "128";

    // process optional arguments
    for (int i = 1; i < argc - requiredArguments; i++) {
        const string& arg = argv[i];

        if (arg == "-s" || arg == "--sa-sparseness") {
            if (i + 1 < argc) {
                saSparse = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-c" || arg == "--cp-sparseness") {
            if (i + 1 < argc) {
                cpSparse = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        }
    }

    string baseFile = argv[argc - 2];
    uint k = atoi(argv[argc - 1]);
    length_t saSF = stoi(saSparse);
    if (saSF == 0 || saSF > 256 || (saSF & (saSF - 1)) != 0) {
        cerr << saSF
             << " is not allowed as sparse factor, should be in 2^[0, 8]"
             << endl;
        return EXIT_FAILURE;
    }
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

    cout << "Welcome to Nexus Stats!\n";

    FMIndexDBG<FMPos> bwt(baseFile, saSF, cpSF, k, false);

    uint32_t numberOfNodes = 0;
    uint64_t numberOfEdges = 0;
    uint64_t totalLength = 0;
    uint64_t totalMultiplicity = 0;
    vector<uint32_t> lengths;
    vector<uint32_t> multiplicities;

    bwt.stats(numberOfNodes, numberOfEdges, totalLength, totalMultiplicity,
              lengths, multiplicities);

    uint32_t medianLength;

    sort(lengths.begin(), lengths.end());
    if (lengths.size() % 2 == 0) {
        medianLength =
            (lengths[lengths.size() / 2 - 1] + lengths[lengths.size() / 2]) / 2;
    } else {
        medianLength = lengths[lengths.size() / 2];
    }

    uint32_t medianMultiplicity;

    sort(multiplicities.begin(), multiplicities.end());
    if (multiplicities.size() % 2 == 0) {
        medianMultiplicity = (multiplicities[multiplicities.size() / 2 - 1] +
                              multiplicities[multiplicities.size() / 2]) /
                             2;
    } else {
        medianMultiplicity = multiplicities[multiplicities.size() / 2];
    }

    cout << "Total no. graph nodes: " << numberOfNodes << "\n";
    cout << "Total no. graph edges: " << numberOfEdges << "\n";
    cout << "Total node multiplicity: " << totalMultiplicity << "\n";
    cout << "Average node multiplicity: "
         << totalMultiplicity / (double)(numberOfNodes) << endl;
    cout << "Median node multiplicity: " << medianMultiplicity << endl;
    cout << "Total node length: " << totalLength << "\n";
    cout << "Average node length: " << totalLength / (double)(numberOfNodes)
         << endl;
    cout << "Median node length: " << medianLength << endl;

    cout << "Bye...\n";
}
