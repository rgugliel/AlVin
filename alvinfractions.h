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
 * \file alvinfractions.h
 * \author Rafael Guglielmetti
 *
 * \class AlVinFractions
 * \brief This class represents a set of possible fractions x_0^2 / (e,e)
 * We generate series of the type:
 * 	(x0 + y)^2 / (e_max, e_max), ..., x0^2 / 1
 */

#ifndef ALVINFRACTIONS_H
#define ALVINFRACTIONS_H

#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

#include "algebraicinteger.h"
#include "alvinfraction.h"

class AlVinFractions {
protected:
  vector<AlgebraicInteger *> possibleNorms2; ///< Possible values for (e,e)
  AlgebraicInteger *possibleNorms2Max;

  vector<AlVinFraction *> alvinFractions;
  vector<AlVinFraction *>::const_iterator
      alvinfractions_it; ///< Iterator to the next element to be returned

  unsigned int lastMaximum;
  unsigned int batchSize;

public:
  /*! \fn AlVinFractions(const unsigned int& iPossibleNorm2)
   * 	\brief Normal constructor
   * 	\param iPossibleNorm2(const unsigned int&): Possible values for (e,e)
   */
  AlVinFractions(vector<AlgebraicInteger *> possibleNorms2);

  virtual ~AlVinFractions();

  /*! \fn getNextAlVinFraction()
   * 	\brief Return the next fractions
   * 	\return Vector of fractions (we can have several couple (x_0^2,(e,e))
   * with same quotient)
   */
  vector<AlVinFraction *> getNextAlVinFraction();

  vector<AlgebraicInteger *> get_possibleNorm2() const;
  const vector<AlgebraicInteger *> *get_ptrPossibleNorm2() const;

private:
  virtual void computeNextAlVinFractions() = 0;
};

#endif // ALVINFRACTIONS_H
