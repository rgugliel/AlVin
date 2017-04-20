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
 * \file infinitensymetries.h
 * \author Rafael Guglielmetti
 * 
 * \class InfiniteNSymetries
 * \brief Try to find integral symmetries of the polyhedron which do not have any common fixed point inside the hyperbolic space. If success: the form is not reflective.
*/

#ifndef INFINITENSYMETRIES_H
#define INFINITENSYMETRIES_H

#include "alvin.h"
#include "CoxIter/lib/numbers/rational.h"

#include <eigen3/Eigen/Dense>
#include <igraph/igraph.h>

#ifdef _OPENMP
#include <omp.h>
#else
inline unsigned int omp_get_thread_num() { return 0; }
inline unsigned int omp_get_max_threads() { return 1; }
#endif

using namespace Eigen;

///< \struct GraphInvolution infinitensymetries.h "One involution corresponding to one integral transformation preserving the polyhedron of H^n"
struct GraphInvolution
{
	vector< unsigned int > iVertices; ///< Vertices of definition
	vector< unsigned int > iPermutation; ///< Permutation defining the transformation
	vector< unsigned int > iLinearlyIndependantVertices; ///< Strictly speaking, the indices i (in iVertices!) such that the set vectors[ iVertices[i] ] is linearly independant
	unsigned int iOrbitCount; ///< Number of invariant subsets (in iVertices)
};

class InfiniteNSymetries
{
	protected:
		AlVin* alvin;
		vector< AlgebraicInteger* > aiQF;
		
		vector< vector< unsigned int > > iGraphMatrix; ///< 1 if bold, 2 if dotted, weight otherwise
		vector< vector< unsigned int > > iCoxeterMatrix;
		
		unsigned int iVectorsCount; ///< Number of vectors computed
		unsigned int iDimension; ///< Dimension of the hyperbolic space
		const unsigned int iVectorSize; ///< iDimension + 1
		
		unsigned int iFixedPointsDimension; ///< Actual dimension of the space of fixed points
		bool bFinished; ///< If true, the form is non-reflective
		
		vector< GraphInvolution > usefulInvolutions;
		
	public:
		InfiniteNSymetries( AlVin* alvin );
		bool Run( const unsigned int& iNRMin, const unsigned int& iNRMax );
		
		unsigned int get_iFixedPointsDimension() const;
		virtual void print_basisFixedPoints( const string& strSpacer = "" ) const = 0;
		
		virtual ~InfiniteNSymetries();
		
		vector< GraphInvolution > get_usefulInvolutions() const;
		
	protected:
		virtual bool bDottedSameWeight( const unsigned int& v1, const unsigned int& v2, const unsigned int& w1, const unsigned int& w2 ) const = 0;
		virtual bool FindIntegralSymmetryFromSubgraph( const vector< unsigned int >& iVertices ) = 0;
};

#endif // INFINITENSYMETRIES_H
