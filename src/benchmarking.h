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

#ifndef BENCHMARKING_H
#define BENCHMARKING_H

#include "searchstrategy.h"

using namespace std;

extern vector<string> schemes;

template <class T>
void printMatches(vector<TextOccurrence> matches, string text, bool printLine,
                  string duration, FMIndex<T>& mapper, string name);

string getFileExt(const string& s);

vector<pair<string, string>> getReads(const string& file);

double
avgVec(vector<length_t> const& v); // note: the average must not be an integer

length_t sum(vector<length_t> const& v);

void writeToOutputSFI(
    const string& file,
    const vector<std::map<std::vector<uint32_t>,
                          std::vector<TextOccurrenceSFI>>>& mPerRead,
    const vector<pair<string, string>>& reads);

void writeToOutputSFR(const string& file,
                      const std::vector<std::vector<FMOcc<FMPosSFR>>>& mPerRead,
                      const vector<pair<string, string>>& reads,
                      FMIndexDBG<FMPosSFR>& index);

double findMedian(vector<length_t> a, int n);

template <class T, class positionClass>
double doBenchSFI(vector<pair<string, string>>& reads, T& mapper,
                  SearchStrategyDBG<T, positionClass>* strategy,
                  string readsFile, length_t ED, std::string outputFile = "");

double doBenchSFR(vector<pair<string, string>>& reads,
                  FMIndexDBG<FMPosSFR>& mapper,
                  SearchStrategyDBG<FMIndexDBG<FMPosSFR>, FMPosSFR>* strategy,
                  string readsFile, length_t ED, std::string cpSparse,
                  std::string outputFile = "");

#endif