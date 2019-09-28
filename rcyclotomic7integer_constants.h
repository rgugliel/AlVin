/*
Copyright (C) 2014, 2015, 2016
Rafael Guglielmetti, rafael.guglielmetti@unifr.ch
*/

/*
This file is part of AlVin.

CoxIter is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

CoxIter is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with AlVin. If not, see <http://www.gnu.org/licenses/>.
*/

/*!
 * \file rcyclotomic7integer_constants.h
 * \author Rafael Guglielmetti
 *
 * \brief Some decompositions of rational primes into RC7
 */

#ifndef RCYCLOTOMIC7INTEGER_CONSTANTS_H
#define RCYCLOTOMIC7INTEGER_CONSTANTS_H

#include "rcyclotomic7integer.h"

#include <algorithm>
#include <gaol/gaol.h>
#include <map>

class RCyclotomic7Integer;

namespace RCyclotomic7Integer_constants {
extern const double dl1; ///< 2 * cos( 2 * pi / 7 )
extern const double dl2; ///< 2 * cos( 4 * pi / 7 )
extern const double dl3; ///< 2 * cos( 6 * pi / 7 )

extern const interval gaol_l1; ///< 2 * cos( 2 * pi / 7 )
extern const interval gaol_l2; ///< 2 * cos( 4 * pi / 7 )
extern const interval gaol_l3; ///< 2 * cos( 6 * pi / 7 )

extern map<unsigned int, vector<RCyclotomic7Integer>>
    iPrimesDecomposition; ///< A table of decomposition of some primes
                          ///< (typically primes up to 10000)
extern const unsigned int iPrimesDecomposition_max;
} // namespace RCyclotomic7Integer_constants

#endif // RCYCLOTOMIC7INTEGER_CONSTANTS_H
