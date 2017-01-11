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
 * \file alvinfractions.h
 * \author Rafael Guglielmetti
 * 
 * \class AlVinFractions
 * \brief This class represents a set of possible fractions x_0^2 / (e,e)
 * We generate series of the type:
 * 	(x0 + y)^2 / (e_max, e_max), ..., x0^2 / 1
*/


#ifndef ALVINFRACTIONS_H
#define ALVINFRACTIONS_H

#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

#include "algebraicinteger.h"
#include "alvinfraction.h"

class AlVinFractions
{
	protected:
		vector< AlgebraicInteger* > aiPossibleNorms2; ///< Possible values for (e,e)
		AlgebraicInteger* aiPossibleNorms2_max;
		
		vector< AlVinFraction* > alvinfractions;
		vector< AlVinFraction* >::const_iterator alvinfractions_it; ///< Iterator to the next element to be returned
		
		unsigned int iLastMaximum;
		unsigned int iBatchSize;
		
	public:
		/*! \fn AlVinFractions( const unsigned int& iPossibleNorm2_max )
		 * 	\brief Normal constructor
		 * 	\param iPossibleNorm2_max( const unsigned int& ): Biggest value for (e,e)
		 */
		AlVinFractions( vector< AlgebraicInteger* > aiPossibleNorms2 );
		
		virtual ~AlVinFractions();
		
		/*! \fn getNextAlVinFraction()
		 * 	\brief Return the next fractions
		 * 	\return Vector of fractions (we can have several couple (x_0^2,(e,e)) with same quotient)
		 */
		vector< AlVinFraction* > getNextAlVinFraction();
		
		vector< AlgebraicInteger* > get_aiPossibleNorm2() const;
		const vector< AlgebraicInteger* >* get_ptraiPossibleNorm2() const;
		
	private:
		virtual void computeNextAlVinFractions() = 0;
};

#endif // ALVINFRACTIONS_H
