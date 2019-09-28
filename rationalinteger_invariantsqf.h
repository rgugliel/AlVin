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
 * \file rationalinteger_invariantsqf.h
 * \author Rafael Guglielmetti
 *
 * \class InvariantsQF
 * \brief Computation of the commensurability of a quadratic form over the
 * rationals
 */

#ifndef INVARIANTSQF_H
#define INVARIANTSQF_H

/*!
 * \file rationalinteger_invariantsqf.h
 * \author Rafael Guglielmetti
 *
 * \class InvariantsQF
 * \brief Compute the invariants of an integral quadratic form of signature
 * (n,1)
 */

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "CoxIter/lib/math_tools.h"
#include "CoxIter/lib/string.h"

using namespace std;
using namespace MathTools;

class InvariantsQF {
private:
  vector<int> iQF;
  unsigned int iDimension;

  vector<vector<unsigned int>> iQF_primes;

  vector<unsigned int> iRamification;
  int iSignedDeterminant;
  vector<unsigned int> iPrimesDeterminant;
  int iDeterminant;
  vector<unsigned int> iRamifiedPrimesQuadraticExtension;

  string strInvariant;

public:
  InvariantsQF(const vector<int> &iQuadraticFormCoeffs);

  vector<unsigned int> get_iRamification() const;
  string get_strInvariant() const;

private:
  void computeInvariants();

  void primeNumbers();
  void computeRamification();

  /*!	\func iRamificationQuaternionAlgebra
   * 	Find the ramification of an elementary quaternion algebra with
   * coefficients in Z
   *
   * 	\param a( int ) -1 or a prime number or -prime
   * 	\param b( int ) -1 or a prime number or -prime
   */
  vector<unsigned int> iRamificationElementaryQuaternionAlgebra(int a,
                                                                int b) const;

  vector<unsigned int> iRamificationProduct(const vector<unsigned int> &iR1,
                                            const vector<unsigned int> &iR2);
};

#endif // INVARIANTSQF_H
