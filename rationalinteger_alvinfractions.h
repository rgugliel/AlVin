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
 * \file rationalinteger_alvinfractions.h
 * \author Rafael Guglielmetti
 *
 * \class RationalInteger_VFs
 * \brief Enumerations of fractions
 */

#ifndef RATIONAL_INTEGER_ALVINFRACTIONS_H
#define RATIONAL_INTEGER_ALVINFRACTIONS_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

#include "alvinfractions.h"
#include "rationalinteger.h"

class RationalInteger_VFs : public AlVinFractions {
private:
  unsigned iPossibleNorms2_max;

public:
  RationalInteger_VFs(vector<AlgebraicInteger *> aiPossibleNorms2);
  ~RationalInteger_VFs();

private:
  void computeNextAlVinFractions();
};

#endif // RATIONAL_INTEGER_ALVINFRACTIONS_H
