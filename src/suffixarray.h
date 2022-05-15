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

#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H

#include "bitvec.h"
#include <fstream> // I/O

#include <string>
#include <vector>

typedef uint64_t length_t;

class SparseSuffixArray {
  private:
    length_t sparseNessFactor;
    Bitvec bitvector;
    std::vector<length_t> sparseSA;

  public:
    bool operator[](const length_t i) const {
        return bitvector[i];
    }

    length_t getSparsenessFactor() const {
        return sparseNessFactor;
    }

    SparseSuffixArray() : sparseNessFactor(0) {
    }

    length_t get(const length_t i) const {
        assert(bitvector[i]);
        return sparseSA[bitvector.rank(i)];
    }

    SparseSuffixArray(const std::vector<length_t>& sa,
                      const length_t sparseNessFactor)
        : sparseNessFactor(sparseNessFactor) {
        bitvector = Bitvec(sa.size());
        sparseSA.reserve(sa.size() / sparseNessFactor);
        for (length_t i = 0; i < sa.size(); i++) {
            const auto& el = sa[i];
            if (el % sparseNessFactor == 0) {
                sparseSA.emplace_back(el);
                bitvector[i] = true;
            }
        }
        bitvector.index();
    }

    SparseSuffixArray(const std::string& basename, const length_t sparseNess)
        : sparseNessFactor(sparseNess) {
        {

            auto name = basename + ".sa.bv." + std::to_string(sparseNessFactor);
            std::ifstream ifs(name);
            if (!ifs) {
                throw std::runtime_error("Cannot open file: " + name +
                                         "\nPick a sparseness factor for which "
                                         "the suffix array was built.");
            }
            bitvector.read(ifs);
        }
        {
            auto name = basename + ".sa." + std::to_string(sparseNessFactor);
            std::ifstream ifs(name);
            if (!ifs) {
                throw std::runtime_error("Cannot open file: " + name +
                                         "\nPick a sparseness factor for which "
                                         "the suffix array was built.");
            }
            ifs.seekg(0, std::ios::end);
            sparseSA.resize(ifs.tellg() / sizeof(length_t));
            ifs.seekg(0, std::ios::beg);
            ifs.read((char*)&sparseSA[0], sparseSA.size() * sizeof(length_t));
        }
    }

    void write(const std::string& basename) const {
        {
            std::ofstream ofs(basename + ".sa.bv." +
                              std::to_string(sparseNessFactor));
            bitvector.write(ofs);
        }
        {
            std::ofstream ofs(basename + ".sa." +
                              std::to_string(sparseNessFactor));
            ofs.write((char*)sparseSA.data(),
                      sparseSA.size() * sizeof(length_t));
        }
    }
};

#endif