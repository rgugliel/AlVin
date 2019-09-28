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
 * \file rationalinteger.h
 * \author Rafael Guglielmetti
 *
 * \class RationalInteger
 * \brief Rational integers
 */

#ifndef RATIONALINTEGER_H
#define RATIONALINTEGER_H

#include "algebraicinteger.h"

class RationalInteger : public AlgebraicInteger // TODO: tout passer en long
{
public:
  long int iVal; // TODO: version avec gmp?

public:
  RationalInteger();
  RationalInteger(long int iVal);
  RationalInteger(const RationalInteger &ri);
  virtual ~RationalInteger();

  AlgebraicInteger *copy() const;
  AlgebraicInteger *aiCopyToInteger(const int &n) const;
  virtual void set(const int &n);
  virtual void set(AlgebraicInteger *ai);

  virtual void removeSquareFactors();

  virtual bool bIsInvertible() const;
  virtual bool bIsSquareOfIvertible() const;

  virtual void gcd(const AlgebraicInteger *ai);
  virtual bool bIsDivisbleBy(const AlgebraicInteger *) const;
  virtual bool divideByIfDivisible(const AlgebraicInteger *ai);
  virtual void divideBy(const AlgebraicInteger *ai);
  virtual void multiplyBy(const int &n);
  virtual void multiplyBy(const AlgebraicInteger *ai);
  virtual void add(const AlgebraicInteger *ai);
  virtual void substract(const AlgebraicInteger *ai);
  virtual void opp();

  virtual bool bIsLessThan(const int &n) const;
  virtual bool bIsLessThan(const AlgebraicInteger &ai) const;
  virtual bool bIsLessOEThan(const AlgebraicInteger &ai) const;
  virtual bool bIsGreaterThan(const int &n) const;
  virtual bool bIsGreaterOEThan(const int &n) const;
  virtual bool bIsEqualTo(const AlgebraicInteger &ai) const;
  virtual bool bIsEqualTo(const int &n) const;

  ostream &print(ostream &) const;
  virtual string to_string(const string &strFormat = "generic",
                           const bool &bProtect = false) const;
  virtual double to_double() const;
  virtual string get_classname() const;

  long int get_iValue() const;

  // --------------------------------------------
  // Operators
  RationalInteger &operator=(const RationalInteger &);
  RationalInteger &operator/=(RationalInteger const &ri);
  RationalInteger &operator*=(RationalInteger const &ri);
  RationalInteger operator-() const;

  RationalInteger operator+(const RationalInteger &) const;
  RationalInteger operator*(const RationalInteger &)const;
  RationalInteger operator-(const RationalInteger &) const;

  bool operator==(const RationalInteger &) const;
  bool operator==(const long int &) const;
  bool operator>(const RationalInteger &) const;
  bool operator<(const RationalInteger &) const;
};

#endif // RATIONALINTEGER_H
