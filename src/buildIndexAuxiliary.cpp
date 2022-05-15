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

#include "buildIndexAuxiliary.h"
#include "bwtrepr.h"
#include "suffixarray.h"

using namespace std;

void readText(const string& filename, string& buf) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    buf.resize(ifs.tellg());
    ifs.seekg(0, ios::beg);
    ifs.read((char*)buf.data(), buf.size());
}

void readSATextMode(const string& filename, vector<length_t>& sa,
                    size_t saSizeHint) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    sa.reserve(saSizeHint);
    length_t el;
    while (ifs >> el)
        sa.push_back(el);
}

void readSA(const string& filename, vector<length_t>& sa, size_t saSizeHint) {
    ifstream ifs(filename, ios::binary);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    size_t numElements = ifs.tellg() / sizeof(length_t);

    if (numElements == saSizeHint) { // file is likely binary
        sa.resize(ifs.tellg() / sizeof(length_t));
        ifs.seekg(0, ios::beg);
        ifs.read((char*)sa.data(), sa.size() * sizeof(length_t));
    } else { // try to read SA in text mode
        readSATextMode(filename, sa, saSizeHint);
    }
}

void sanityCheck(const string& T, vector<length_t>& sa) {
    // check T for correctness
    if (T.back() == '\n')
        throw runtime_error("T should end with a \'$\' character, "
                            "not with a newline");

    if (T.back() != '$')
        throw runtime_error("T should end with a \'$\' character");

    if (sa.size() != T.size())
        throw runtime_error("Text and suffix array contain a "
                            "different number of elements");

    // briefly check the suffix array
    length_t min = *min_element(sa.begin(), sa.end());
    length_t max = *max_element(sa.begin(), sa.end());

    if (min == 1 && max == T.size()) { // rebase to [0..T.size()-1]
        for (auto& el : sa)
            el--;
        min--;
        max--;
    }

    if (min != 0 || max != T.size() - 1)
        throw runtime_error("Suffix array must contain numbers between "
                            "[0 and " +
                            to_string(T.size() - 1) + "]");

    // check if all numbers in the suffix array are present
    Bitvec bv(sa.size());
    for (length_t i : sa)
        bv[i] = true;

    for (size_t i = 0; i < bv.size(); i++)
        if (!bv[i])
            throw runtime_error("Suffix " + to_string(i) +
                                " seems "
                                "to be missing from suffix array");

    // extra check:
    //      we could check T to see if the SA correctly sorts suffixes of T
}
