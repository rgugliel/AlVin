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
 * \file alvinfraction.h
 * \author Rafael Guglielmetti
 * 
 * \class AlVinFraction
 * \brief This class represents one fraction x0^2 / (e,e)
*/

#ifndef ALVINFRACTION_H
#define ALVINFRACTION_H

#include "algebraicinteger.h"

class AlVinFraction
{
	public: // For simplicity, variables are public be are meant to be read-only. Careful.
		AlgebraicInteger* aiX0; ///< x_0
		AlgebraicInteger* aiNorm2; ///< (e,e)
		AlgebraicInteger* aiNumerator; ///< When normalized with respect to the biggest possible value for (e,e): x_0^2 / (e,e) = ( x_0^2 * (e_max, e_max)/(e,e) ) / (e_max, e_max)
		
	public:
		AlVinFraction( AlgebraicInteger* aiX0, AlgebraicInteger* aiNorm2, AlgebraicInteger* aiNumerator );
		~AlVinFraction();
		
		bool operator==( AlVinFraction const& ) const;
		bool operator!=( AlVinFraction const& ) const;
		
		friend bool operator<( const AlVinFraction &f1, const AlVinFraction &f2 );
		friend ostream& operator<<( ostream& , AlVinFraction const & );
};

/*!	\fn isLessThanPtrAlVinFraction( const AlVinFraction* a, const AlVinFraction* b )
 * 	\brief Check if *a < *b
*/
bool isLessThanPtrAlVinFraction( const AlVinFraction* a, const AlVinFraction* b );

#endif // ALVINFRACTION_H
