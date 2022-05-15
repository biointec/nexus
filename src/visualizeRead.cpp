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
    cout << "Usage: ./visualizeRead [options] basefilename read\n\n";
    cout << " [Pattern matching options]\n";
    cout << " -sfr --strain-free\tstrain-free matching\n";
    cout << " -e   --max-ed\t\tmaximum edit distance [default = 0]\n";
    cout << " -s   --sa-sparseness\tsuffix array sparseness factor "
            "[default = "
            "1]\n";
    cout << " -c   --cp-sparseness\tsparseness factor that indicates "
            "how many checkpoints must be stored to identify nodes. Use "
            "\"none\" to use no checkpoints. Choose a value that was also used "
            "during the building process. "
            "[default = 128]\n";
    cout
        << " -f   --filter\t\tfiltering type that should be used to filter the "
           "occurrences. This option is only valid in case of strain-free "
           "matching. Options:\n\t"
        << "linear\t\tlinear filtering is efficient but does not filter out "
           "all redundant occurrences. Additionally, in some exceptional "
           "cases, a non-optimal replacement occurrence can be chosen. This "
           "is the default option.\n\t"
        << "complete\tcomplete filtering leads to a set of occurrences with "
           "no redundancy. This option is very slow however and thus not "
           "recommended.\n";
    cout << " -p   --partitioning\t\tAdd flag to do uniform/static/dynamic "
            "partitioning. Dynamic partitioning cannot be used with "
            "strain-free matching. [default = static]\n";
    cout << " -m   --metric\t\tAdd flag to set distance metric "
            "(editnaive/editopt/hamming) [default = editopt]\n";
    cout << " -ss  --search-scheme\tChoose the search scheme\n  options:\n\t"
         << "kuch1\tKucherov k + 1\n\t"
         << "kuch2\tKucherov k + 2\n\t"
         << "kianfar\tOptimal Kianfar scheme\n\t"
         << "manbest\tManual best improvement for Kianfar scheme (only for ed "
            "= 4)\n\t"
         << "pigeon\tPigeonhole scheme\n\t"
         << "01*0\t01*0 search scheme\n\t"
         << "naive\tnaive backtracking\n\t"
         << "custom\tcustom search scheme, the next parameter should be a path "
            "to the folder containing this search scheme\n\n";

    cout << " [Visualization options]\n";
    cout << " -d   --visualization-depth\t\tDepth of the visualized "
            "neighborhood around the paths of interest [default = 3]\n";
    cout << " -o   --output-files\t\tPrefix of the output files that will be "
            "created during the visualization process [default = "
            "basefilename]\n\n";

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

    string saSparse = "1";
    string cpSparse = "128";
    string maxED = "0";
    string visDepthString = "3";
    string searchscheme = "kuch1";
    string customFile = "";
    string outputFile = "";
    bool strainFree = false;
    bool filteringIsChosen = false;
    bool filteringOptionComplete = false;

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
        } else if (arg == "-d" || arg == "--visualization-depth") {
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
            return false;
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

    string baseFile = argv[argc - 2];
    string read = argv[argc - 1];

    if (outputFile == "") {
        outputFile = baseFile;
    }

    cout << "Start creation of BWT approximate matcher on graphs" << endl;

    if (strainFree) {

        FMIndexDBG<FMPosSFR> bwt(baseFile, saSF, cpSF, strainFree,
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
        bwt.visualizeSubgraphs(results, visDepth, outputFile);

        delete strategy;

    } else {

        FMIndexDBG<FMPos> bwt(baseFile, saSF, cpSF, strainFree);

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
        bwt.visualizeSubgraphs(results, visDepth, outputFile);

        delete strategy;
    }

    cout << "Bye...\n";
}
