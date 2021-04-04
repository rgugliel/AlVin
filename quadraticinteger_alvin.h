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
 * \file quadraticinteger_alvin.h
 * \author Rafael Guglielmetti
 *
 * \class QuadraticInteger_AlVin
 * \brief AlVin for quadratic integers
 */

#ifndef QUADRATICINTEGER_ALVIN_H
#define QUADRATICINTEGER_ALVIN_H

#include "alvin.h"
#include "quadraticinteger.h"
#include "quadraticinteger_alvinfractions.h"

class QuadraticInteger_AlVin : public AlVin {
private:
  vector<QuadraticInteger> qiQF;  ///< Local copy of the coefficients of the
                                  ///< quadratic form; used to work faster
  vector<QuadraticInteger> qi2QF; ///< Coefficients of the quadratic form * 2
  vector<QuadraticInteger> qiDQF; ///< Coefficients of the quadratic form * d
  vector<QuadraticInteger>
      qiDQFConj; ///< Conjugates of coefficients of the quadratic form * d

  vector<QuadraticInteger>
      qiBilinearProducts; ///< To control during findVector that the vector has
                          ///< negative product with the preceding ones
  vector<QuadraticInteger> qiVectorCurrent;     ///< Current working vector
  vector<vector<QuadraticInteger *>> qiVectors; ///< List of vectors found

  long int d;      ///< Local copy of QuadraticInteger::d
  long int dISqrt; ///< Local copy of isqrt( QuadraticInteger::d )
  bool bIsOneMod4; ///< Local copy of QuadraticInteger::bIsOneMod4

public:
  QuadraticInteger_AlVin(const vector<QuadraticInteger> &iQuadraticFormCoeffs,
                         const string &strOuputMathematicalFormat,
                         const bool &bWriteInfo, const bool &bDebug);

  virtual string get_strField() const;

private:
  void findPossibleNorms2();

  bool PreRun();

  void findVector(AlgebraicInteger *aiX0, AlgebraicInteger *aiNorm2);

  /*!	\fn findVector( QuadraticInteger* qi0, QuadraticInteger* qiNorm2,
   * unsigned int iIndex, QuadraticInteger qiSumComp, QuadraticInteger
   * qiGCDComponents ); \brief Find a potential vector The vector e must satisfy
   * the cristallographic condition, (e,e) = iNorm2 and the first component is
   * specified Remark: we do not run other tests
   *
   * 	\param qi0( QuadraticInteger* ): First component of the vector
   * 	\param qiNorm2( QuadraticInteger* ): Square of the norm of the vector
   * 	\param iIndex( const unsigned int& ): Current component of the vector
   * 	\param qiSumComp( const unsigned int& ): iSumComp =
   * \sum_{i=iIndex}^iDimension a_i * x_i^2 \param qiGCDComponents(
   * QuadraticInteger ): GCD of the found components
   */
  void findVector(QuadraticInteger *qi0, QuadraticInteger *qiNorm2,
                  unsigned int iIndex, QuadraticInteger qiSumComp,
                  QuadraticInteger qiGCDComponents);

  void addCandidate();
  virtual void addVectorChild(const vector<AlgebraicInteger *> &aiVector);
  virtual int addVector_findWeight(AlgebraicInteger *aiNumerator,
                                   AlgebraicInteger *aiDenominator);
};

#endif // QUADRATICINTEGER_ALVIN_H
