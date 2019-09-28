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
 * \file rationalinteger_alvin.h
 * \author Rafael Guglielmetti
 *
 * \class RationalInteger_AlVin
 * \brief Find the vectors for rational integers
 */

#ifndef RATIONAL_INTEGER_ALVIN_H
#define RATIONAL_INTEGER_ALVIN_H

#include <iterator>

#include "alvin.h"
#include "rationalinteger.h"
#include "rationalinteger_alvinfractions.h"

class RationalInteger_AlVin : public AlVin {
private:
  vector<unsigned int> riQF;  ///< Local copy of the coefficients of the
                              ///< quadratic form; used to work faster
  vector<unsigned int> iGCDs; ///< We have [ i ] = gcd( a_i, a_i+1, ..., a_n )
                              ///< or 0 if the gcd equals to one
  bool bLastTwoCoefficientsQFAre1Mod4; ///< True if the last two coefficients of
                                       ///< the quadratic form are 1 mod 4

  vector<int>
      iBilinearProducts; ///< To control during findVector that the vector has
                         ///< negative product with the preceding ones
  vector<int> iVectorCurrent;   ///< Current working vector
  vector<vector<int>> iVectors; ///< List of vectors found

public:
  RationalInteger_AlVin(const vector<int> &iQuadraticFormCoeffs,
                        const string &strOuputMathematicalFormat,
                        const bool &bWriteInfo, const bool &bDebug);

  virtual string get_strField() const;

private:
  void findPossibleNorms2();

  bool PreRun();

  void findVector(AlgebraicInteger *aiX0, AlgebraicInteger *aiNorm2);

  /*!	\fn findVector( const unsigned int& i0, const unsigned int& iNorm2,
   * unsigned int iIndex = 1, unsigned int iSumComp = 0, unsigned int
   * iGCDComponents = 1 ); \brief Find a potential vector The vector e must
   * satisfy the cristallographic condition, (e,e) = iNorm2 and the first
   * component is specified Remark: we do not run other tests
   *
   * 	\param i0( const unsigned int& ): First component of the vector
   * 	\param iNorm2( const unsigned int& ): Square of the norm of the vector
   * 	\param iIndex( const unsigned int& ): Current component of the vector
   * 	\param iSumComp( const unsigned int& ): iSumComp =
   * \sum_{i=iIndex}^iDimension a_i * x_i^2 \param iGCDComponents( unsigned int
   * ): GCD of the found components of the vector
   */
  void findVector(const unsigned int &i0, const unsigned int &iNorm2,
                  unsigned int iIndex = 1, unsigned int iSumComp = 0,
                  unsigned int iGCDComponents = 1);

  void addCandidate();
  bool bCandidatePreserveEvenLattice(const unsigned int &iNorm2) const;

  virtual void addVectorChild(const vector<AlgebraicInteger *> &aiVector);
  virtual int addVector_iFindWeight(AlgebraicInteger *aiNumerator,
                                    AlgebraicInteger *aiDenominator);
};

#endif // RATIONAL_INTEGER_ALVIN_H
