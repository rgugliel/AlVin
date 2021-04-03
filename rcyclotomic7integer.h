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
 * \file rcyclotomic7integer.h
 * \author Rafael Guglielmetti
 *
 * \class RCyclotomic7Integer
 * \brief RC7 and their operations
 */

#ifndef RCYCLOTOMI7CINTEGER_H
#define RCYCLOTOMI7CINTEGER_H

#include <algorithm>
#include <gaol/gaol.h>
#include <gmpxx.h>
#include <map>
#include <pari/pari.h>
#include <string>

#include "algebraicinteger.h"
#include "rcyclotomic7integer_constants.h"

class RCyclotomic7Integer : public AlgebraicInteger {
public:
  array<mpz_class, 3>
      iC; ///< Decomposition in the basis lambda_1, ..., lambda_(iDegree)

private:
  static bool bClassInitialized; ///< If the class has been initialized

  static vector<vector<int>>
      iPowersL; ///< Decomposition of l^2, l^3, ..., l^iDegree
  static array<array<array<mpq_class, 2>, 3>, 2>
      mpqRationalApproximations; ///< Rational approximations (10^-50, 10^-200)

public:
  RCyclotomic7Integer(const RCyclotomic7Integer &rci);
  RCyclotomic7Integer(const int &n);
  RCyclotomic7Integer(const mpz_class &n);
  RCyclotomic7Integer(const array<mpz_class, 3> &iCoefficients);
  RCyclotomic7Integer();
  virtual ~RCyclotomic7Integer();

  virtual void set(AlgebraicInteger *ai);
  virtual void set(const int &n);
  AlgebraicInteger *copy() const;
  AlgebraicInteger *aiCopyToInteger(const int &n) const;

  virtual void removeSquareFactors();

  virtual bool isSquareOfIvertible() const;
  virtual bool isInvertible() const;
  bool bIsAssociateTo(const AlgebraicInteger *ai);

  virtual void gcd(const AlgebraicInteger *ai);
  virtual bool isDivisbleBy(const AlgebraicInteger *ai) const;
  virtual void divideBy(const AlgebraicInteger *ai);
  virtual bool divideByIfDivisible(const AlgebraicInteger *ai);
  virtual void multiplyBy(const AlgebraicInteger *ai);
  virtual void multiplyBy(const int &n);
  virtual void opp();
  virtual void substract(const AlgebraicInteger *ai);
  virtual void add(const AlgebraicInteger *ai);

  virtual bool isEqualTo(const int &n) const;
  virtual bool isEqualTo(const AlgebraicInteger &ai) const;
  virtual bool isGreaterOEThan(const int &n) const;
  virtual bool isGreaterThan(const int &n) const;
  virtual bool isLessOEThan(const AlgebraicInteger &ai) const;
  virtual bool isLessThan(const AlgebraicInteger &ai) const;
  virtual bool isLessThan(const int &n) const;

  ostream &print(ostream &) const;
  virtual std::string get_classname() const;

  virtual std::string to_string(const string &strFormat = "generic",
                                const bool &bProtect = false) const;
  virtual double to_double() const;
  interval to_interval() const;

  long int floor() const;

public:
  static void initialize();
  mpz_class iNorm() const;

  bool bIsLongInt() const;

  void conjugate(unsigned int i);

  vector<RCyclotomic7Integer> rciPrimeFactors() const;
  map<RCyclotomic7Integer, unsigned int> rciPrimeDecomposition() const;

  static void initializePrimeDecomposition();
  static bool computePrimeDecomposition(const unsigned int &iPrime);

  // --------------------------------------------
  // Operators
  RCyclotomic7Integer &operator=(const RCyclotomic7Integer &rci);
  RCyclotomic7Integer &operator/=(RCyclotomic7Integer const &rci);
  RCyclotomic7Integer &operator*=(RCyclotomic7Integer const &rci);

  RCyclotomic7Integer operator+(const RCyclotomic7Integer &) const;
  RCyclotomic7Integer operator*(const RCyclotomic7Integer &)const;
  RCyclotomic7Integer operator-(const RCyclotomic7Integer &) const;
  RCyclotomic7Integer operator-() const;

  bool operator>(const RCyclotomic7Integer &rci) const;
  bool operator==(const RCyclotomic7Integer &rci) const;
  bool operator==(const long int &) const;

private:
  static void computePrimesDecomposition_factorMinimalPolynomial(
      GEN &gFactors, const unsigned int &iPrime);
  static void
  computePrimesDecomposition_addFactor(const unsigned int &iPrime,
                                       const unsigned int &iFactorsCount,
                                       RCyclotomic7Integer rci);
};

bool operator<(const RCyclotomic7Integer &rci1,
               const RCyclotomic7Integer &rci2);

inline interval gaol_intervalFromMPZclass(const mpz_class &i) {
  if (i.fits_slong_p())
    return interval(i.get_si());
  else
    return interval(i.get_str().c_str());
}

inline interval gaol_SQRTQuotient(RCyclotomic7Integer rciNum,
                                  const RCyclotomic7Integer &rciDenom) {
  if (rciNum.divideByIfDivisible(&rciDenom)) {
    if (rciNum.bIsLongInt())
      return sqrt(interval(-rciNum.iC[0].get_si()));
    else
      return sqrt(rciNum.to_interval());
  } else
    return sqrt(rciNum.to_interval() / rciDenom.to_interval());
}

#endif // RCYCLOTOMIC7INTEGER_H
