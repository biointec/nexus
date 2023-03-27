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
#include "toACGT.h"

#include <dirent.h>
#include <sstream>

using namespace std;

typedef uint64_t length_t;

void showUsage() {
    cout
        << "This program constructs a new implicit pan-genome de Bruijn\n"
           "graph, along with the underlying bidirectional FM-index.\n\n\n"

           "Usage: ./nexusBuild <base filename> <k_list>\n\n"

           " Following input parameters are required:\n"
           "  <base filename>       base filename of the output index and\n"
           "                        possibly the input text if it was already\n"
           "                        preprocessed\n"
           "  <k_list>              the de Bruijn parameter: a "
           "comma-separated\n"
           "                        list of integers is required (e.g.,\n"
           "                        20,21,23)\n\n"

           " Following input files are required:\n"
           "  <base filename>.txt:  Nexus v1.1.0 or lower requires a "
           "preprocessed\n"
           "                        input text with all genomes readily "
           "concatenated,\n"
           "                        containing only the following characters:\n"
           "                        A, C, G, T, % and $ (at the very end). No\n"
           "                        newlines are allowed. Nexus v1.1.1 or "
           "higher\n"
           "                        can build the index from .fasta files "
           "directly.\n\n\n"

           " [options]\n"
           "  --skip                Skip the building process of the data\n"
           "                        structures that are independent of the de\n"
           "                        Bruijn k parameter (i.e., the "
           "bidirectional\n"
           "                        FM-index). These data structures must be\n"
           "                        available in the directory.\n\n"
           "  -f/--fasta            Build the index directly from\n"
           "                        .fasta/.fa/.fna files. Contigs or\n"
           "                        chromosomes within the fasta files are\n"
           "                        ignored and concatenated together. The\n"
           "                        input fasta files should be added "
           "directly\n"
           "                        after this option. Multiple files can be\n"
           "                        added. Alternatively, if no filenames are\n"
           "                        passed on, Nexus scans the working\n"
           "                        directory and includes all fasta files it\n"
           "                        can find.\n\n"
           "  --fastaWithC          This option is analogous to the --fasta\n"
           "                        option, but instead of concatenating\n"
           "                        contigs or chromosomes together, they are\n"
           "                        considered to be different strains in the\n"
           "                        pan-genome.\n\n"
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
           "  -p/--progress         Report extra progress updates\n\n\n";
}

/**
 * @brief Check if a string ends with a certain suffix
 *
 * @param str The string to be investigated
 * @param suffix The suffix to be found
 * @return true if the string ends with the suffix
 * @return false otherwise
 */
bool endsWith(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() &&
           0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

/**
 * @brief Check if a string starts with a certain prefix
 *
 * @param str The string to be investigated
 * @param prefix The prefix to be found
 * @return true if the string starts with the prefix
 * @return false otherwise
 */
bool startsWith(const std::string& str, const std::string& prefix) {
    return str.size() >= prefix.size() &&
           0 == str.compare(0, prefix.size(), prefix);
}

/**
 * @brief Add an additional strain from a fasta file to the pan-genome
 *
 * @param fastaFile the fasta file to extract the strain from
 * @param referenceText the reference text or pan-genome to add the strain to
 * @param retainContigs bool indicating wether the contigs/chromosomes in the
 * fasta file should be retained as strains or concatenated
 */
void updateReferenceText(string& fastaFile, string& referenceText,
                         const bool& retainContigs) {
    cout << "Preprocessing " << fastaFile << "... ";
    ifstream ifs(fastaFile);
    if (!ifs)
        throw runtime_error("Cannot open file: " + fastaFile);
    std::string line;
    bool firstLine = true;
    while (std::getline(ifs, line)) {
        if (!(line[0] == '>')) {
            // Keep only A, C, G or T characters on the line
            toACGT(line);
            referenceText += line;
        } else if (retainContigs && !firstLine) {
            // If contigs should be considered as different strains, then add a
            // separation character when a header line is encountered that was
            // not the first one
            referenceText += "%";
        } else {
            // Keep track of whether we are on the first header line or not
            firstLine = false;
        }
    }
    ifs.close();
    // Add a separation character to mark the end of this fasta file
    referenceText += "%";
    cout << "Done." << endl;
}

bool parseArguments(int argc, char* argv[], string& baseFN, vector<uint>& k,
                    vector<int>& saSF, vector<int>& cpSF, bool& progress,
                    bool& skip, string& referenceText) {
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
        } else if (((arg == "-f" || (arg == "--fasta") ||
                     (arg == "--fastaWithC")))) {
            bool retainContigs = false;
            if (arg == "--fastaWithC") {
                retainContigs = true;
            }
            i++;
            arg = (argv[i]);
            // Loop over subsequent arguments as long as they point to fasta
            // files
            while (endsWith(arg, ".fasta") || endsWith(arg, ".fa") ||
                   endsWith(arg, ".fna")) {
                // Add the strain in the file to the reference text
                updateReferenceText(arg, referenceText, retainContigs);
                // Go to the next argument
                i++;
                arg = (argv[i]);
            }
            // Check if the reference text is empty, which means that no fasta
            // file arguments were passed specifically
            if (referenceText.empty()) {
                // In this case, we loop over the working directory to find all
                // fasta files and include them in the pan-genome.
                DIR* dir;
                struct dirent* ent;
                if ((dir = opendir(".")) != NULL) {
                    /* find all the files and directories within directory */
                    while ((ent = readdir(dir)) != NULL) {
                        string fastaFile = ent->d_name;
                        // Check if it is a fasta file
                        if (endsWith(fastaFile, ".fasta") ||
                            endsWith(fastaFile, ".fa") ||
                            endsWith(fastaFile, ".fna")) {
                            // Add the strain(s) to the pan-genome
                            updateReferenceText(fastaFile, referenceText,
                                                retainContigs);
                        }
                    }
                    closedir(dir);
                } else {
                    /* could not open directory */
                    throw runtime_error("Could not open the working directory");
                }
            }
            // Check if we effectively found a pan-genome
            if (referenceText.empty()) {
                throw runtime_error(
                    "No fasta files were found for the reference text");
            }
            // Replace the last % separation character with the sentinel
            // character
            referenceText[referenceText.size() - 1] = '$';
            // Reset the i iterator such that we do not skip an option
            i--;
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
    string referenceText = "";

    if (!parseArguments(argc, argv, baseFilename, k, saSF, cpSF, progress, skip,
                        referenceText)) {
        showUsage();
        return EXIT_FAILURE;
    }

    if (!referenceText.empty()) {
        // Print the number of occurrences of each nucleotide in the pan-genome
        printCharacterSubstitutionState();
    }

    cout << "Welcome to Nexus!\n";
    cout << "Alphabet size is " << ALPHABET - 1 << " + 1\n";
    // cout << "k is " << k << "\n";

    try {
        // cout << "Start creation of BWT approximate matcher" << endl;

        FMIndexDBG<FMPos>::buildFMIndexDBG(baseFilename, k, saSF, cpSF,
                                           progress, skip, referenceText);
    } catch (const std::exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    cout << "Exiting... bye!" << endl;
    return EXIT_SUCCESS;
}