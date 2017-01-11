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
 * \file notreflexive.h
 * \author Rafael Guglielmetti
 * 
 * \class NotReflexive
 * \brief Create systems of equations to test the non-reflexivity of a quadratic form defined over Z
*/

#ifndef NOTREFLEXIVE_H
#define NOTREFLEXIVE_H

#include "alvin.h"

///< \struct NotReflexive_Graph notreflexive.h "Structure which corresponds to one Euclidean graph which cannot be extended"
struct NotReflexive_Graph
{
	vector< short unsigned int > iGraphVertices; ///< Vertices of the graph
	vector< short unsigned int > iVariablesName; ///< For each vertex, the corresponding variable
	vector< short unsigned int > iVariablesToCoeff; ///< [ i ] = j means that x_i has the jth coefficient in the quadratic form
	
	vector< short unsigned int > iVariablesGreaterThan; ///< [ i ] = j means x_i >= x_j
	vector< AlgebraicInteger* > aiVariablesCount; ///< Number of time each variable appear (useful for the norm equation)
};

class NotReflexive
{
	protected:
		AlVin *alvin;
		vector< vector< AlgebraicInteger* > > aiVectors;
		unsigned int iDimension;
		
		vector< AlgebraicInteger* > aiQF;
		vector< AlgebraicInteger* > ai2QF;
		
		vector< vector< NotReflexive_Graph > > graphs; ///< The first index is for the number of variables, then one for each graph which cannot be extended
		
		string strOFormat;
		string strAlgebraicIntegerType;
		
		vector< AlgebraicInteger* > aiPossibleNorm2;
		
	public:
		NotReflexive( AlVin *v );
		virtual ~NotReflexive();
		
		void Run();
		
	private:
		void prepareGraphsList();
		void createSystemsEquations();
		virtual void createSystemEquations( NotReflexive_Graph nrg ) = 0;
};

#endif // NOTREFLEXIVE_H
