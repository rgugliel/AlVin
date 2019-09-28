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
 * \file rcyclotomic7integer_alvin.h
 * \author Rafael Guglielmetti
 *
 * \class RCyclotomic7Integer_AlVin
 * \brief AlVin for RC7
 */

#ifndef RCYCLOTOMIC7INTEGER_ALVIN_H
#define RCYCLOTOMIC7INTEGER_ALVIN_H

#include "alvin.h"
#include "rcyclotomic7integer.h"
#include "rcyclotomic7integer_alvinfractions.h"

class RCyclotomic7Integer_AlVin : public AlVin {
private:
  vector<RCyclotomic7Integer> rciQF; ///< Local copy of the coefficients of the
                                     ///< quadratic form; used to work faster
  vector<RCyclotomic7Integer>
      rci2QF; ///< Coefficients of the quadratic form * 2

  vector<RCyclotomic7Integer>
      rciBilinearProducts; ///< To control during findVector that the vector has
                           ///< negative product with the preceding ones
  vector<RCyclotomic7Integer> rciVectorCurrent;     ///< Current working vector
  vector<vector<RCyclotomic7Integer *>> rciVectors; ///< List of vectors found

public:
  RCyclotomic7Integer_AlVin(
      const vector<RCyclotomic7Integer> &rciQuadraticFormCoeffs,
      const string &strOuputMathematicalFormat, const bool &bWriteInfo,
      const bool &bDebug);

  virtual std::string get_strField() const;
  virtual bool PreRun();

private:
  void findPossibleNorms2();
  void findVector(AlgebraicInteger *aiX0, AlgebraicInteger *aiNorm2);
  void addVectorChild(const vector<AlgebraicInteger *> &aiVector);
  virtual int addVector_iFindWeight(AlgebraicInteger *aiNumerator,
                                    AlgebraicInteger *aiDenominator);

  void findVector(RCyclotomic7Integer *rci0, RCyclotomic7Integer *rciNorm2,
                  unsigned int iIndex, RCyclotomic7Integer rciSumComp,
                  RCyclotomic7Integer rciGCDComponents);
  void addCandidate();

  virtual void print_initialInformationChild() const;
};

#endif // RCYCLOTOMIC7INTEGER_ALVIN_H
