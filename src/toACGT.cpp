#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <set>
#include <vector>

using namespace std;

// a map containing the counts for the characters in the preprocessed pan-genome
map<char, int64_t> newCharacterCounts = {
    {'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};

// a map containing the counts for the characters in the original strains
map<char, int64_t> originalCharacterCounts = {
    {'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}, {'R', 0},
    {'Y', 0}, {'M', 0}, {'K', 0}, {'S', 0}, {'W', 0},
    {'H', 0}, {'B', 0}, {'V', 0}, {'D', 0}, {'N', 0}};

// random generator
static std::mt19937 gen(4);

/**
 * @brief Select a random element from a set
 *
 * @param start Start of the set
 * @param end End of the set
 * @return Iter - selected element
 */
template <typename Iter> Iter select_randomly(Iter start, Iter end) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(gen));
    return start;
}

void printCharacterSubstitutionState() {

    std::cout << "\nOriginal character counts in the pan-genome strains: "
              << endl;
    for (auto it = originalCharacterCounts.cbegin();
         it != originalCharacterCounts.cend(); ++it) {
        if (it->second) {
            std::cout << it->first << ": " << it->second << "\n";
        }
    }

    std::cout << "\nNew character counts in the pan-genome strains: " << endl;
    for (auto it = newCharacterCounts.cbegin(); it != newCharacterCounts.cend();
         ++it) {
        std::cout << it->first << ": " << it->second << "\n";
    }
    std::cout << endl;
}

void toACGT(string& line) {
    // Substitution sets

    // puRine
    std::vector<char> R = {'A', 'G'};
    // pYrimidine
    std::vector<char> Y = {'C', 'T'};
    // aMino group on base
    std::vector<char> M = {'A', 'C'};
    // Keto group on base
    std::vector<char> K = {'G', 'T'};
    // Strong base pairing
    std::vector<char> S = {'C', 'G'};
    // Weak base pairing
    std::vector<char> W = {'A', 'T'};
    // Not G
    std::vector<char> H = {'A', 'C', 'T'};
    // Not A
    std::vector<char> B = {'C', 'G', 'T'};
    // Not T
    std::vector<char> V = {'A', 'C', 'G'};
    // Not C
    std::vector<char> D = {'A', 'G', 'T'};
    // aNy base
    std::vector<char> N = {'A', 'C', 'G', 'T'};

    char newChar;

    for (size_t i = 0; i < line.size(); i++) {
        // Capitalize
        char ch = toupper(line[i]);
        originalCharacterCounts[ch]++;
        switch (ch) {
        case 'A':
            newCharacterCounts['A']++;
            break;
        case 'C':
            newCharacterCounts['C']++;
            break;
        case 'G':
            newCharacterCounts['G']++;
            break;
        case 'T':
            newCharacterCounts['T']++;
            break;
        case 'R':
            newChar = *select_randomly(R.begin(), R.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'Y':
            newChar = *select_randomly(Y.begin(), Y.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'M':
            newChar = *select_randomly(M.begin(), M.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'K':
            newChar = *select_randomly(K.begin(), K.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'S':
            newChar = *select_randomly(S.begin(), S.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'W':
            newChar = *select_randomly(W.begin(), W.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'H':
            newChar = *select_randomly(H.begin(), H.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'B':
            newChar = *select_randomly(B.begin(), B.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'V':
            newChar = *select_randomly(V.begin(), V.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'D':
            newChar = *select_randomly(D.begin(), D.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case 'N':
            newChar = *select_randomly(N.begin(), N.end());
            line[i] = newChar;
            newCharacterCounts[newChar]++;
            break;
        case '\n':
            break;
        default:
            string errorMessage = "Unknown character: no replacement can be "
                                  "chosen for the following character: ";
            errorMessage.push_back(ch);
            throw runtime_error(errorMessage);
        }
    }
}