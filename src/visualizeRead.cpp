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

vector<string> schemes = {"kuch1",  "kuch2", "kianfar", "manbest",
                          "pigeon", "01*0",  "custom",  "naive"};

/**
 * @brief Show the usage in terminal
 *
 */
void showUsage() {
    cout <<

        "This program aligns a short, single end read to a pan-genome de\n"
        "Bruijn graph. It visualizes the corresponding node paths with their\n"
        "surroundings in the graph.\n\n\n"

        "Usage: ./visualizeRead [options] <basefilename> <k> <read>\n\n"

        " Following input parameters are required:\n"
        "  <basefilename>      base filename of the input index\n"
        "  <k>                 the de Bruijn parameter of the index\n"
        "  <read>              the read that must be aligned and "
        "visualized.\n\n\n"

        " [options]\n"
        "  -e/--max-ed         maximum edit distance [default = 0]\n\n"

        "  -s/--sa-sparseness  suffix array sparseness factor [default = "
        "16]\n\n"

        "  -c/--cp-sparseness  sparseness factor that indicates how many\n"
        "                      checkpoints must be stored to identify nodes.\n"
        "                      Use \"none\" to use no checkpoints. Choose a\n"
        "                      value that was also used during the building\n"
        "                      process. [default = 128]\n\n"

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
        "searchscheme\n\n"

        "  -sfr/--strain-free  strain-free matching: occurrences can be\n"
        "                      identified as any path of connected nodes. In\n"
        "                      other words, they do not have to occur exactly\n"
        "                      in one of the input genomes of the pan-genome.\n"
        "                      This is option is not activated by default and\n"
        "                      is slower than the default implementation.\n\n"

        "  -f/--filter         filtering type that should be used to filter\n"
        "                      the occurrences. This option is only valid in\n"
        "                      case of strain-free matching. Options:\n"
        "                       * linear: linear filtering is efficient but\n"
        "                         does not filter out all redundant\n"
        "                         occurrences. Additionally, in some\n"
        "                         exceptional cases, a non-optimal "
        "replacement\n"
        "                         occurrence can be chosen. This is the\n"
        "                         default option.\n"
        "                       * complete: complete filtering leads to a set\n"
        "                         of occurrences with no redundancy. This\n"
        "                         option is very slow however and thus not\n"
        "                         recommended.\n\n"

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

    int requiredArguments = 3; // baseFile of files, k and file containing reads

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

    string saSparse = "16";
    string cpSparse = "128";
    string maxED = "0";
    string visDepthString = "3";
    string searchscheme = "kuch1";
    string customFile = "";
    string outputFile = "";
    bool strainFree = false;
    bool filteringIsChosen = false;
    bool filteringOptionComplete = false;
    bool separateEdges = true;

    PartitionStrategy pStrat = STATIC;
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
                    if (strainFree == true) {
                        throw runtime_error(
                            "Dynamic partitioning cannot be used with "
                            "strain-free matching.");
                    }
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
        } else if (arg == "-f" || arg == "--filter") {
            filteringIsChosen = true;
            if (i + 1 < argc) {
                string s = argv[++i];
                if (s == "linear") {
                    filteringOptionComplete = false;
                } else if (s == "complete") {
                    filteringOptionComplete = true;
                } else {
                    throw runtime_error(
                        s + " is not a filtering option\nOptions are: "
                            "linear, complete");
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
        } else if (arg == "-d" || arg == "--depth") {
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
                    if (strainFree == true) {
                        throw runtime_error(
                            "Hamming distance is currently not an option in "
                            "combination with strain-free matching.");
                    }
                } else {
                    throw runtime_error(s +
                                        " is not a metric option\nOptions are: "
                                        "editopt, editnaive, hamming");
                }

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-sfr" || arg == "--strain-free") {
            strainFree = true;
            if (metric == HAMMING) {
                throw runtime_error(
                    "Hamming distance is currently not an option in "
                    "combination with strain-free matching.");
            }
            if (pStrat == DYNAMIC) {
                throw runtime_error("Dynamic partitioning cannot be used with "
                                    "strain-free matching.");
            }
        }

        else {
            cerr << "Unknown argument: " << arg << " is not an option" << endl;
            return EXIT_FAILURE;
        }
    }

    if (filteringIsChosen && !strainFree) {
        throw runtime_error("The filtering option should not be set if "
                            "strain-free matching is not chosen.");
    }

    length_t ed = stoi(maxED);
    if (ed < 0 || ed > 5) {
        cerr << ed << " is not allowed as maxED should be in [0, 5]" << endl;

        return EXIT_FAILURE;
    }
    length_t saSF = stoi(saSparse);
    if (saSF == 0 || saSF > 256 || (saSF & (saSF - 1)) != 0) {
        cerr << saSF
             << " is not allowed as sparse factor, should be in 2^[0, 8]"
             << endl;
        return EXIT_FAILURE;
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

    if (ed != 4 && searchscheme == "manbest") {
        throw runtime_error("manbest only supports 4 allowed errors");
    }

    string baseFile = argv[argc - 3];
    uint k = atoi(argv[argc - 2]);
    string read = argv[argc - 1];

    if (outputFile == "") {
        outputFile = baseFile;
    }

    cout << "Start creation of BWT approximate matcher on graphs" << endl;

    if (strainFree) {

        FMIndexDBG<FMPosSFR> bwt(baseFile, saSF, cpSF, k, strainFree,
                                 filteringOptionComplete);

        SearchStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>* strategy;
        if (searchscheme == "kuch1") {
            strategy = new KucherovKplus1DBG<FMIndexDBG<FMPosSFR>, FMPosSFR>(
                bwt, pStrat, metric);
        } else if (searchscheme == "kuch2") {
            strategy = new KucherovKplus2DBG<FMIndexDBG<FMPosSFR>, FMPosSFR>(
                bwt, pStrat, metric);
        } else if (searchscheme == "kianfar") {
            strategy = new OptimalKianfarDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>(
                bwt, pStrat, metric);
        } else if (searchscheme == "manbest") {
            strategy = new ManBestStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>(
                bwt, pStrat, metric);
        } else if (searchscheme == "01*0") {
            strategy =
                new O1StarSearchStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>(
                    bwt, pStrat, metric);
        } else if (searchscheme == "pigeon") {
            strategy =
                new PigeonHoleSearchStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>(
                    bwt, pStrat, metric);
        } else if (searchscheme == "custom") {
            strategy =
                new CustomSearchStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>(
                    bwt, customFile, pStrat, metric);
        } else if (searchscheme == "naive") {
            strategy =
                new NaiveBackTrackingStrategyDBG<FMIndexDBG<FMPosSFR>,
                                                 FMPosSFR>(bwt, pStrat, metric);
        } else {
            // should not get here
            throw runtime_error(searchscheme +
                                " is not on option as search scheme");
        }
        StrainFreeMapper mapper(strategy);
        auto results = mapper.matchApproxSFR(read, ed);
        bwt.visualizeSubgraphs(results, visDepth, outputFile, separateEdges);

        delete strategy;

    } else {

        FMIndexDBG<FMPos> bwt(baseFile, saSF, cpSF, k, strainFree);

        SearchStrategyDBG<FMIndexDBG<FMPos>, FMPos>* strategy;
        if (searchscheme == "kuch1") {
            strategy = new KucherovKplus1DBG<FMIndexDBG<FMPos>, FMPos>(
                bwt, pStrat, metric);
        } else if (searchscheme == "kuch2") {
            strategy = new KucherovKplus2DBG<FMIndexDBG<FMPos>, FMPos>(
                bwt, pStrat, metric);
        } else if (searchscheme == "kianfar") {
            strategy = new OptimalKianfarDBG<FMIndexDBG<FMPos>, FMPos>(
                bwt, pStrat, metric);
        } else if (searchscheme == "manbest") {
            strategy = new ManBestStrategyDBG<FMIndexDBG<FMPos>, FMPos>(
                bwt, pStrat, metric);
        } else if (searchscheme == "01*0") {
            strategy = new O1StarSearchStrategyDBG<FMIndexDBG<FMPos>, FMPos>(
                bwt, pStrat, metric);
        } else if (searchscheme == "pigeon") {
            strategy =
                new PigeonHoleSearchStrategyDBG<FMIndexDBG<FMPos>, FMPos>(
                    bwt, pStrat, metric);
        } else if (searchscheme == "custom") {
            strategy = new CustomSearchStrategyDBG<FMIndexDBG<FMPos>, FMPos>(
                bwt, customFile, pStrat, metric);
        } else if (searchscheme == "naive") {
            strategy =
                new NaiveBackTrackingStrategyDBG<FMIndexDBG<FMPos>, FMPos>(
                    bwt, pStrat, metric);
        } else {
            // should not get here
            throw runtime_error(searchscheme +
                                " is not on option as search scheme");
        }

        auto results = strategy->matchApproxSFI(read, ed);
        bwt.visualizeSubgraphs(results, visDepth, outputFile, separateEdges);

        delete strategy;
    }

    cout << "Bye...\n";
}
