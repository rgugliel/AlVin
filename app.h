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
 * \file app.h
 * \author Rafael Guglielmetti
 * 
 * \class App
 * \brief Main class
*/

#ifndef APP_H
#define APP_H

#include "rationalinteger_alvin.h"
#include "rationalinteger_notreflexive.h"
#include "rationalinteger_invariantsqf.h"
#include "rationalinteger_infinitensymetries.h"

#include "quadraticinteger_alvin.h"
#include "quadraticinteger_infinitensymetries.h"

#ifdef _RC7AVAILABLE_
	#include "rcyclotomic7integer_alvin.h"
	#include "rcyclotomic7integer_infinitensymetries.h"
#endif

#include <iostream>
#include <regex>
#include <chrono>

using namespace std;

class App
{
	private:
		vector< AlgebraicInteger* > aiQF; ///< Coefficients of the quadratic form
		
		bool bCheckNR; ///< If true, try to determine if the form is non reflexive (needed: maxv, nrmin, nrmax)
		unsigned int iNRMin; ///< Minimal number of vertices for the subsets
		unsigned int iNRMax; ///< Maximal number of vertices for the subsets
		
		bool bCheckNREquations; ///< Write the system of equations to check if the form is non reflexive
		
		bool bDebug; ///< If true, more information are displayed
		bool bComputeInvariantsQF; ///< If true, compute the invariants of the qf
		bool bComputeInvariantsPolyhedron; ///< If true, compute the invariants of the polyhedra
		int iCreateImage; ///< -1: not specified (i.e. yes if possible and if the number of vectors is <= 25), 0: no, 1: force yes
				
		unsigned int iMinVectors; ///< Minimal number of vectors to compute
		unsigned int iMaxVectors; ///< Maximal number of vectors to compute
		
		string strField; ///< Field of definition (rationals, quadratic)
		unsigned int iFieldSupp; ///< Additional data (for example d if k=Q[ sqrt d ])
		
		string strOuputMathematicalFormat; ///< Format for mathematical output (generic, mathematica)
		
	public:
		App();
		~App();
		
		void readMainParameters( int argc, char **argv );
		void Run();
		
	private:
		AlVin* instanciateAlVin();
		NotReflexive* instanciateNotReflexiveEquations( AlVin* v );
		InfiniteNSymetries* instanciateInfiniteNSymetries( AlVin* v );
};

#endif // APP_H
