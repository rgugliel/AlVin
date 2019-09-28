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
 * \file rcyclotomic7integer_infinitensymetries.h
 * \author Rafael Guglielmetti
 *
 * \class RCyclotomic7Integer_InfiniteNSymetries
 * \brief To find integral symmetries of the space
 */

#ifndef RCYCLOTOMIC7INTEGER_INFINITENSYMETRIES_H
#define RCYCLOTOMIC7INTEGER_INFINITENSYMETRIES_H

#include "infinitensymetries.h"
#include "rcyclotomic7integer.h"

namespace Eigen {
template <class> struct NumTraits;
template <>
struct NumTraits<Rational<RCyclotomic7Integer>> ///< \struct
                                                ///< NumTraits<Rational<RCyclotomic7Integer>>
                                                ///< rcyclotomic7integer_infinitensymetries.h
                                                ///< "Specialization for Eigen
                                                ///< matrices to
                                                ///< RCyclotomic7Integer type"
{
  typedef Rational<RCyclotomic7Integer> Real;
  typedef Rational<RCyclotomic7Integer> NonInteger;
  typedef Rational<RCyclotomic7Integer> Nested;
  typedef Rational<RCyclotomic7Integer> Literal;

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
struct significant_decimals_impl<
    Rational<RCyclotomic7Integer>> ///< \struct
                                   ///< significant_decimals_impl<Rational<RCyclotomic7Integer>>
                                   ///< rcyclotomic7integer_infinitensymetries.h
                                   ///< "Specialization for Eigen matrices to
                                   ///< RCyclotomic7Integer type"
{
  // Infinite precision when printing
  static inline int run() { return 0; }
};
} // namespace internal
}; // namespace Eigen

class RCyclotomic7Integer_InfiniteNSymetries : public InfiniteNSymetries {
private:
  vector<Matrix<Rational<RCyclotomic7Integer>, Dynamic, 1>>
      rciVectorsC; ///< Vectors found by AlVin, in column
  vector<vector<Rational<RCyclotomic7Integer>>>
      qiDottedWeights; ///< Weights of the dotted lines

  Matrix<Rational<RCyclotomic7Integer>, Dynamic, Dynamic>
      rciBasisFixedPoints; ///< Basis of the fixed points of the involutions, in
                           ///< column

  vector<Rational<RCyclotomic7Integer>>
      rciQF; ///< Coefficients of the diagonal quadratic form (in fractions)
  vector<vector<RCyclotomic7Integer>>
      riVectorsProducts;               ///< Products of each pair of vecteurs
  vector<long int> iVectorsNormsFloor; ///< Length of vectors

public:
  RCyclotomic7Integer_InfiniteNSymetries(AlVin *alvin);
  ~RCyclotomic7Integer_InfiniteNSymetries();

  Rational<RCyclotomic7Integer> rciVectorNorm(
      const Matrix<Rational<RCyclotomic7Integer>, Dynamic, 1> &rciV) const;
  Rational<RCyclotomic7Integer> rciVectorsProduct(
      const Matrix<Rational<RCyclotomic7Integer>, Dynamic, 1> &rciV1,
      const Matrix<Rational<RCyclotomic7Integer>, Dynamic, 1> &rciV2) const;
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
