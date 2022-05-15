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

#include "fmocc.h"

// Definition of function pointers
replacesPtr FMOccSFR::replaces = &FMOccSFR::replacesRegular;
setReplacementPtr FMOccSFR::setReplacement = &FMOccSFR::setReplacementRegular;
reverseNodePathPtr FMOccSFR::reverseNodePath =
    &FMOccSFR::reverseNodePathRegular;
reportPtr FMOccSFR::report = &FMOccSFR::reportRegular;

// ============================================================================
// CLASS BIFMOCC
// ============================================================================

std::ostream& operator<<(std::ostream& o, const FMOcc<FMPos>& m) {
    return o << m.getPosition() << "\tEdit distance: " << m.getDistance()
             << "\tShift: " << m.getShift();
}

std::ostream& operator<<(std::ostream& o, const FMOcc<FMPosSFR>& m) {
    return o << m.getPosition() << "\tEdit distance: " << m.getDistance()
             << "\tShift: " << m.getShift();
}