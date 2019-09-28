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
 * \file rcyclotomic7integer_alvinfractions.h
 * \author Rafael Guglielmetti
 *
 * \class RCyclotomic7Integer_VFs
 * \brief Enumerations of fractions
 */

#ifndef RCYCLOTOMIC7INTEGER_VFS_H
#define RCYCLOTOMIC7INTEGER_VFS_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

#include "alvinfractions.h"
#include "rcyclotomic7integer.h"

class RCyclotomic7Integer_VFs : public AlVinFractions {
private:
  RCyclotomic7Integer rciAlpha0; ///< First coefficient of the quadratic form
  RCyclotomic7Integer
      rciMinusAlpha0; ///< -First coefficient of the quadratic form
  RCyclotomic7Integer *rciPossibleNorms2_max;

public:
  RCyclotomic7Integer_VFs(vector<AlgebraicInteger *> aiPossibleNorms2,
                          const RCyclotomic7Integer &rciAlpha0);
  ~RCyclotomic7Integer_VFs();

private:
  void computeNextAlVinFractions();
};

#endif // RCYCLOTOMIC7INTEGER_VFS_H
