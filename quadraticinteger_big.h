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
 * \file quadraticinteger_big.h
 * \author Rafael Guglielmetti
 *
 * \class QuadraticIntegerBig
 * \brief Quadratic integers with bigint components
 */

#ifndef QUADRATICINTEGER_BIG_H
#define QUADRATICINTEGER_BIG_H

#include "algebraicinteger.h"
#include "quadraticinteger.h"

#include <array>
#include <cmath>
#include <gmpxx.h>
#include <map>
#include <string>
#include <vector>

using namespace std;

class QuadraticIntegerBig : public AlgebraicInteger {
public:
  mpz_class
      a; ///< First component of the quadratic integer (in the usual Z-basis)
  mpz_class
      b; ///< Second component of the quadratic integer (in the usual Z-basis)

  static int d;             ///< We work in Q[ sqrt d ]
  static mpf_class sqrtd;   ///< Square root of d
  static bool bIsOneMod4;   ///< True if d is 1 mod 4, false otherwise
  static int iDiscriminant; ///< Discriminant of the quadratic field

  static map<unsigned int, vector<unsigned int>>
      iPellMinimalSolution; ///< Solutions for the Pell equations (used to
                            ///< compute the decompositions of rational prime
                            ///< numbers)
  static map<unsigned int, vector<long int>>
      iFundamentalUnits; ///< Fundamental units

public:
  QuadraticIntegerBig();
  QuadraticIntegerBig(const int &iVal);
  QuadraticIntegerBig(const mpz_class &iVal);
  QuadraticIntegerBig(const QuadraticIntegerBig &qi);
  QuadraticIntegerBig(const QuadraticInteger &qi);
  QuadraticIntegerBig(const int &a, const int &b);
  QuadraticIntegerBig(const mpz_class &a, const mpz_class &b);
  virtual ~QuadraticIntegerBig();

  AlgebraicInteger *copy() const;
  AlgebraicInteger *aiCopyToInteger(const int &n) const;
  virtual void set(const int &n);
  virtual void set(AlgebraicInteger *ai);

  virtual void removeSquareFactors();

  virtual bool isInvertible() const;
  virtual bool isSquareOfIvertible() const;

  virtual void gcd(const AlgebraicInteger *ai);
  virtual bool isDivisbleBy(const AlgebraicInteger *) const;
  virtual void divideBy(const AlgebraicInteger *ai);
  virtual bool divideByIfDivisible(const AlgebraicInteger *ai);
  virtual void multiplyBy(const int &n);
  virtual void multiplyBy(const AlgebraicInteger *ai);
  virtual void add(const AlgebraicInteger *ai);
  virtual void substract(const AlgebraicInteger *ai);
  virtual void opp();

  virtual bool isLessThan(const int &n) const;
  virtual bool bIsLessThan(const long int &n) const;
  virtual bool isLessThan(const AlgebraicInteger &ai) const;
  virtual bool isLessOEThan(const AlgebraicInteger &ai) const;
  virtual bool isGreaterThan(const int &n) const;
  virtual bool bIsGreaterThan(const long int &n) const;
  virtual bool isGreaterOEThan(const int &n) const;
  virtual bool isEqualTo(const AlgebraicInteger &ai) const;
  virtual bool isEqualTo(const int &n) const;

  ostream &print(ostream &) const;
  virtual string to_string(const string &strFormat = "generic",
                           const bool &bProtect = false) const;
  virtual double to_double() const;
  virtual string get_classname() const;

  mpz_class floor() const;

public:
  static bool bIsDAdmissible(const unsigned int &d);
  static void set_d(const unsigned int &d);
  static vector<QuadraticIntegerBig>
  qiFactorsRationalPrime(const unsigned int &iPrime,
                         bool bWithMultiplicities = false);
  static array<long int, 2> iPellEquation(const unsigned int &iPrime);
  static mpz_class iSQRT_quotient(const QuadraticIntegerBig &qiNum,
                                  const QuadraticIntegerBig &qiDen);
  static mpz_class iSQRTsup_quotient(const QuadraticIntegerBig &qiNum,
                                     const QuadraticIntegerBig &qiDen);

  void conjugate();
  mpz_class iNorm() const;
  mpz_class iTrace() const;

  int iValuation(const QuadraticIntegerBig &qi);

  vector<QuadraticIntegerBig> qiPrimeFactors() const;
  map<QuadraticIntegerBig, unsigned int> qiPrimeDecomposition() const;

  bool bIsAssociateTo(QuadraticIntegerBig qi2);

  // --------------------------------------------
  // Operators
  QuadraticIntegerBig &operator=(const QuadraticIntegerBig &);
  QuadraticIntegerBig &operator/=(QuadraticIntegerBig const &qi);
  QuadraticIntegerBig &operator*=(QuadraticIntegerBig const &qi);

  QuadraticIntegerBig operator+(const QuadraticIntegerBig &) const;
  QuadraticIntegerBig operator*(const QuadraticIntegerBig &)const;
  QuadraticIntegerBig operator-(const QuadraticIntegerBig &) const;
  QuadraticIntegerBig operator-() const;

  bool operator==(const long int &) const;
  bool operator==(const QuadraticIntegerBig &) const;
  bool operator>(const QuadraticIntegerBig &) const;
};

bool operator<(const QuadraticIntegerBig &qi1, const QuadraticIntegerBig &qi2);

#endif // QUADRATICINTEGER_H
