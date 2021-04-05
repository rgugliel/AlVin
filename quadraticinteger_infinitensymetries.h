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
 * \file quadraticinteger_infinitensymetries.h
 * \author Rafael Guglielmetti
 *
 * \class QuadraticInteger_InfiniteNSymetries
 * \brief To find integral symmetries of the space
 */

#ifndef QUADRATIC_INFINITENSYMETRIES_H
#define QUADRATIC_INFINITENSYMETRIES_H

#include "infinitensymetries.h"
#include "quadraticinteger.h"
#include "quadraticinteger_big.h"

namespace Eigen {
template <class> struct NumTraits;
template <>
struct NumTraits<
    Rational<QuadraticIntegerBig>> ///< \struct
                                   ///< NumTraits<Rational<QuadraticIntegerBig>>
                                   ///< quadraticinteger_infinitensymetries.h
                                   ///< "Specialization for Eigen
                                   ///< matrices to
                                   ///< QuadraticIntegerBig type"
{
  typedef Rational<QuadraticIntegerBig> Real;
  typedef Rational<QuadraticIntegerBig> NonInteger;
  typedef Rational<QuadraticIntegerBig> Nested;
  typedef Rational<QuadraticIntegerBig> Literal;

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
    QuadraticIntegerBig>> ///< \struct
                          ///< significant_decimals_impl<Rational<QuadraticIntegerBig>>
                          ///< quadraticinteger_infinitensymetries.h
                          ///< "Specialization for Eigen matrices to
                          ///< QuadraticIntegerBig type"
{
  // Infinite precision when printing
  static inline int run() { return 0; }
};
} // namespace internal
}; // namespace Eigen

class QuadraticInteger_InfiniteNSymetries : public InfiniteNSymetries {
private:
  vector<Matrix<Rational<QuadraticIntegerBig>, Dynamic, 1>>
      rqiVectorsC; ///< Vectors found by AlVin, in column
  vector<vector<Rational<QuadraticIntegerBig>>>
      qiDottedWeights; ///< Weights of the dotted lines

  Matrix<Rational<QuadraticIntegerBig>, Dynamic, Dynamic>
      rqiBasisFixedPoints; ///< Basis of the fixed points of the involutions, in
                           ///< column

  vector<Rational<QuadraticIntegerBig>>
      rqiQF; ///< Coefficients of the diagonal quadratic form (in fractions)
  vector<vector<QuadraticIntegerBig>>
      riVectorsProducts;              ///< Products of each pair of vecteurs
  vector<long int> vectorsNormsFloor; ///< Length of vectors

public:
  QuadraticInteger_InfiniteNSymetries(AlVin *alvin);
  ~QuadraticInteger_InfiniteNSymetries();

  Rational<QuadraticIntegerBig> rqiVectorNorm(
      const Matrix<Rational<QuadraticIntegerBig>, Dynamic, 1> &rqiV) const;
  Rational<QuadraticIntegerBig> rqiVectorsProduct(
      const Matrix<Rational<QuadraticIntegerBig>, Dynamic, 1> &rqiV1,
      const Matrix<Rational<QuadraticIntegerBig>, Dynamic, 1> &rqiV2) const;
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
