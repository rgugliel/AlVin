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
 * \file alvin.h
 * \author Rafael Guglielmetti
 *
 * \class AlVin
 * \brief Main class for AlVin
 */

#ifndef ALVIN_H
#define ALVIN_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>

#include <map>
#include <memory>

#include "CoxIter/coxiter.h"
#include "CoxIter/lib/math_tools.h"

#include "algebraicinteger.h"
#include "alvinfractions.h"

using namespace std;
using namespace MathTools;

class AlVin {
protected:
  vector<AlgebraicInteger *> aiQF;
  unsigned int iDimension; ///< Dimension of the space

  bool bComputeInvariantsPolyhedron; ///< If yes, we compute the invariants of
                                     ///< the final polyhedron
  bool bDebug;     ///< If true, more information are displayed
  bool bWriteInfo; ///< If we want to write informations (false if AlVin is used
                   ///< "as a plugin")
  int iCreateImage; ///< -1: not specified (i.e. yes if possible and if the
                    ///< number of vectors is <= 25), 0: no, 1: force yes

  string strAlgebraicIntegerType; ///< RationalInteger, QuadraticInteger

  vector<unsigned int> iQBlocksSize; ///< Sizes of blocks of coefficients
  vector<unsigned int>
      iComponentLessThan; ///< The entry [ i ] = j means x_i <= x_j, [ i ] = 0
                          ///< means no information for x_i

  vector<vector<AlgebraicInteger *>> aiVectors; ///< The vectors
  unsigned int iVectorsCount;                   ///< Number of vectors found
  unsigned int iVectorsCount_second; ///< Number of vectors found (second batch)

  vector<vector<AlgebraicInteger *>>
      aiVectors_candidates; ///< Vectors which are compatiable with the previous

  vector<vector<unsigned int>> iCoxeterMatrix;

  vector<AlgebraicInteger *>
      iBilinearProducts; ///< To control during findVector_simple that the
                         ///< vector has negative product with the preceding
                         ///< ones

  AlVinFractions *vf;

  string strOuputMathematicalFormat;

  string strFinalInformation; ///< Final information that will be displayed

private:
  CoxIter *ptrCI;

public:
  /*!	\fn AlVin( const string& strOuputMathematicalFormat, const bool&
   * bWriteInfo, const bool& bDebug ) \brief Constructor for the parent AlVin
   * object \param strOuputMathematicalFormat( const string& ): Format for the
   * output: generic, mathematica, pari, latex \param bWriteInfo( const bool& ):
   * If true, display info \param bDebug( const bool& ): If true, some debug
   * info
   */
  AlVin(const string &strOuputMathematicalFormat, const bool &bWriteInfo,
        const bool &bDebug);

  virtual ~AlVin();

  /*!	\fn bool PreRun() = 0
   * 	\brief Prepare the computations
   */
  virtual bool PreRun() = 0;

  /*!	\fn Run( unsigned int iMinVectors = 0, unsigned int iMaxVectors = 0,
   * bool bLastCheckFV = true ) \brief Run the algo! \param iMinVectors(unsigned
   * int): If specified, does not check the finiteness while this->iVectorsCount
   * < iMinVectors \param iMaxVectors(unsigned int): If specified, stops when
   * this->iVectorsCount == iMaxVectors \param bLastCheckFV(bool): If false,
   * does not check the finiteness for the final polyhedron (useful for debug)
   */
  bool Run(unsigned int iMinVectors = 0, unsigned int iMaxVectors = 0,
           bool bLastCheckFV = true);

  /*!	\fn AlgebraicInteger* aiBilinearProduct( const vector< AlgebraicInteger*
   * >& v1, const vector< AlgebraicInteger* >& v2 ) \brief Compute the product
   * between two vectors
   */
  AlgebraicInteger *aiBilinearProduct(const vector<AlgebraicInteger *> &v1,
                                      const vector<AlgebraicInteger *> &v2);

  /*!	\fn print_iQF
   * 	\brief Print the quadratic form
   */
  void print_iQF() const;

  /*!	\fn print_vectors
   * 	\brief Print the found vectors
   */
  void print_vectors() const;

  /*!	\fn get_iCoxeterMatrix
   * 	\brief Return the Coxeter matrix
   */
  vector<vector<unsigned int>> get_iCoxeterMatrix() const;

  /*!	\fn get_iDimension
   * 	\brief Return n, the dimension
   */
  unsigned int get_iDimension() const;

  /*!	\fn get_strFinalInformation
   * 	\brief Return some information stored during the computations
   */
  string get_strFinalInformation() const;

  /*!	\fn get_strCoxeterMatrix
   * 	\brief Return the Coxeter matrix (string format; rows are separated by a
   * \n)
   */
  string get_strCoxeterMatrix() const;

  /*!	\fn get_strField
   * 	\brief Return the field of definition of the current quadratic form
   */
  virtual string get_strField() const = 0;

  /*!	\fn get_strQF
   * 	\brief Return the current quadratic form (string)
   *
   * 	\param strSeparator( const string& ) Separator between the coefficients
   */
  string get_strQF(const string &strSeparator = ", ") const;

  /*!	\fn string get_strAlgebraicIntegerType() const
   * 	\brief Return the underlying field of the quadratic form
   */
  string get_strAlgebraicIntegerType() const;

  /*!	\fn get_aiVectors
   * 	\brief Return the list of found vectors
   */
  vector<vector<AlgebraicInteger *>> get_aiVectors() const;

  /*!	\fn get_iVectorsCount
   * 	\brief Return the number of found vectors
   */
  unsigned int get_iVectorsCount() const;

  /*!	\fn get_aiQF
   * 	\brief Return the quadratic form
   */
  vector<AlgebraicInteger *> get_aiQF() const;

  /*!	\fn get_aiPossibleNorm2
   * 	\brief Return the list of possible values for (e,e), where e is a root
   * of the quadratic lattice
   */
  vector<AlgebraicInteger *> get_aiPossibleNorm2() const;

  /*!	\fn get_ptraiPossibleNorm2
   * 	\brief Return a pointer to the list of possible values for (e,e), where
   * e is a root of the quadratic lattice
   */
  const vector<AlgebraicInteger *> *get_ptraiPossibleNorm2() const;

  /*!	\fn get_ptrCI
   * 	\brief Return a pointer to the the last CoxIter object
   */
  CoxIter *get_ptrCI() const;

  // ---------------------------------------------------
  // Setters

  /*!	\fn void set_iCreateImage( const int& iValue )
   * 	\brief Create image with the graph: 1: forced yes, 0: forced no, -1:
   * default value (only if at most 25 vectors found)
   */
  void set_iCreateImage(const int &iValue);

  /*!	\fn void set_bComputeInvariantsPolyhedron( const bool& bValue )
   * 	\brief Do we have to compute the invariants of the final polyhedron
   */
  void set_bComputeInvariantsPolyhedron(const bool &bValue);

protected:
  void initializations();

  void print_initialInformation() const;
  virtual void print_initialInformationChild() const;
  void print_finallInformation() const;

private:
  virtual void findPossibleNorms2() = 0;
  virtual void findVector(AlgebraicInteger *aiX0,
                          AlgebraicInteger *aiNorm2) = 0;

  void findFirstVectors();

  void addVector(const vector<AlgebraicInteger *> &aiVect);
  virtual int addVector_iFindWeight(AlgebraicInteger *aiNumerator,
                                    AlgebraicInteger *aiDenominator) = 0;

  void printFoundVector(vector<AlgebraicInteger *> aiV,
                        const unsigned int &iIndex,
                        const bool &bFirst = false) const;

  virtual void addVectorChild(const vector<AlgebraicInteger *> &aiVector) = 0;
};

inline AlgebraicInteger *
AlVin::aiBilinearProduct(const vector<AlgebraicInteger *> &v1,
                         const vector<AlgebraicInteger *> &v2) {
  AlgebraicInteger *aiProduct(aiQF[0]->copy());
  AlgebraicInteger *aiTemp(aiQF[0]->copy());

  aiProduct->opp();
  aiProduct->multiplyBy(v1[0]);
  aiProduct->multiplyBy(v2[0]);

  for (unsigned int j(1); j <= iDimension; j++) {
    aiTemp->set(aiQF[j]);
    aiTemp->multiplyBy(v1[j]);
    aiTemp->multiplyBy(v2[j]);
    aiProduct->add(aiTemp);
  }

  delete aiTemp;

  return aiProduct;
}

#endif // ALVIN_H
