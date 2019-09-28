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
 * \file rationalinteger_infinitensymetries.h
 * \author Rafael Guglielmetti
 *
 * \class RationalInteger_InfiniteNSymetries
 * \brief To find integral symmetries of the space
 */

#ifndef RATIONALINTEGER_INFINITENSYMETRIES_H
#define RATIONALINTEGER_INFINITENSYMETRIES_H

#include "infinitensymetries.h"
#include "rationalinteger.h"

namespace Eigen {
template <class> struct NumTraits;
template <>
struct NumTraits<Rational<
    RationalInteger>> ///< \struct NumTraits<Rational<RationalInteger>>
                      ///< rationalinteger_infinitensymetries.h "Specialization
                      ///< for Eigen matrices to RationalInteger type"
{
  typedef Rational<RationalInteger> Real;
  typedef Rational<RationalInteger> NonInteger;
  typedef Rational<RationalInteger> Nested;
  typedef Rational<RationalInteger> Literal;

  static inline Real epsilon() { return 0; }
  static inline Real dummy_precision() { return 0; }
  enum {
    IsInteger = 0,
    IsSigned = 1,
    IsComplex = 0,
    RequireInitialization = 1,
    ReadCost = Eigen::HugeCost,
    AddCost = Eigen::HugeCost,
    MulCost = Eigen::HugeCost
  };
};

namespace internal {
template <>
struct significant_decimals_impl<Rational<
    RationalInteger>> ///< \struct
                      ///< significant_decimals_impl<Rational<RationalInteger>>
                      ///< rationalinteger_infinitensymetries.h "Specialization
                      ///< for Eigen matrices to RationalInteger type"
{
  // Infinite precision when printing
  static inline int run() { return 0; }
};
} // namespace internal
}; // namespace Eigen

class RationalInteger_InfiniteNSymetries : public InfiniteNSymetries {
private:
  vector<Matrix<Rational<RationalInteger>, Dynamic, 1>>
      rriVectorsC; ///< Vectors found by AlVin, in column
  vector<vector<Rational<RationalInteger>>>
      riDottedWeights; ///< Weights of the dotted lines

  Matrix<Rational<RationalInteger>, Dynamic, Dynamic>
      rriBasisFixedPoints; ///< Basis of the fixed points of the involutions, in
                           ///< column

  vector<Rational<RationalInteger>>
      rriQF; ///< Coefficients of the diagonal quadratic form (in fractions)
  vector<vector<RationalInteger>>
      riVectorsProducts; ///< Products of each pair of vecteurs

public:
  RationalInteger_InfiniteNSymetries(AlVin *alvin);
  ~RationalInteger_InfiniteNSymetries();

  Rational<RationalInteger> rriVectorNorm(
      const Matrix<Rational<RationalInteger>, Dynamic, 1> &rriV) const;
  Rational<RationalInteger> rriVectorsProduct(
      const Matrix<Rational<RationalInteger>, Dynamic, 1> &rriV1,
      const Matrix<Rational<RationalInteger>, Dynamic, 1> &rriV2) const;
  bool FindIntegralSymmetryFromSubgraph(const vector<unsigned int> &iVertices);

  virtual void print_basisFixedPoints(const string &strSpacer = "") const;

private:
  bool bDottedSameWeight(const unsigned int &v1, const unsigned int &v2,
                         const unsigned int &w1, const unsigned int &w2) const;
  void WorkWithIsomorphisms(const vector<unsigned int> &iVertices,
                            const vector<GraphInvolution> &iIsomorphisms);
  vector<GraphInvolution>
  FindIsomorphismsInSubgraph(const vector<unsigned int> &iVertices);
};

#endif // RATIONALINTEGER_INFINITENSYMETRIES_H
