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
 * \file algebraicinteger.h
 * \author Rafael Guglielmetti
 * 
 * \class AlgebraicInteger
 * \brief Parent class for rational, quadratic and rc7 integers
*/

#ifndef ALGEBRAICINTEGER_H
#define ALGEBRAICINTEGER_H

#include <vector>
#include <iostream>
#include <string>

#include "CoxIter/lib/math_tools.h"

using namespace std;
using namespace MathTools;

class AlgebraicInteger
{
	public:
		virtual ~AlgebraicInteger();
		
		/*!	\fn aiCopyToInteger
		 * 	\brief Create an algebraic integer from an integer
		 * 	
		 * 	\param n( const int ): Integer
		 * 	\return Pointer to the algebraic integer created
		 */
		virtual AlgebraicInteger* aiCopyToInteger( const int& n ) const = 0;
		
		/*!	\fn copy
		 * 	\brief Create a copy of an algebraic integer
		 * 	
		 * 	\return Pointer to the algebraic integer created
		 */
		virtual AlgebraicInteger* copy() const = 0;
		
		/*!	\fn set( const int& n ) = 0
		 * 	\brief Assign an integer to the algebraic integer
		 * 	\param n( const int& ) Integer
		 */
		virtual void set( const int& n ) = 0;
		
		/*!	\fn set(  AlgebraicInteger* ai ) = 0
		 * 	\brief Assign an algebraic integer to the algebraic integer
		 * 	\param ai( AlgebraicInteger* ) Pointer to the algebraic integer
		 */
		virtual void set( AlgebraicInteger* ai ) = 0;
		
		/*!	\fn removeSquareFactors()
		 * 	\brief Remove all the square factors of the algebraic integer
		 */
		virtual void removeSquareFactors() = 0;
		
		/*!	\fn bIsInvertible()
		 * 	\brief Check whether the number is invertible
		 * 	\return True if invertible or not
		 */
		virtual bool bIsInvertible() const = 0;
		
		/*!	\fn bIsSquareOfIvertible()
		 * 	\brief Check whether the number is the square of an invertible element
		 * 	\return True if invertible or not
		 */
		virtual bool bIsSquareOfIvertible() const = 0;
		
		/*!	\fn gcd( const AlgebraicInteger* ai )
		 * 	\brief *this becomes the gcd of *this and the parameter
		 * 	\param ai( const AlgebraicInteger* ) The other algebraic integer
		 */
		virtual void gcd( const AlgebraicInteger* ai ) = 0;
		
		/*!	\fn bIsDivisbleBy( const AlgebraicInteger* ai )
		 * 	\brief Check whether *this is divisible by the parameter
		 * 	\param ai( const AlgebraicInteger* ) Pointer to the divisor
		 * 	\return True if *this is divisible by the parameter
		 */
		virtual bool bIsDivisbleBy( const AlgebraicInteger* ai ) const = 0;
		
		/*!	\fn divideBy( const AlgebraicInteger* ai )
		 * 	\brief Performs the division of *this by the parameter
		 * 	Remark: we do not check that the *this is divisible by *ai. If it is not the case, the division will be an approximation!
		 * 	\param ai( const AlgebraicInteger* ) Pointer to the divisor
		 */
		virtual void divideBy( const AlgebraicInteger *ai ) = 0;
		
		/*!	\fn multiplyBy( const int& n )
		 * 	\brief Performs the multiplication of *this by the parameter
		 * 	\param ai( const AlgebraicInteger* ) Pointer to the multiplicand
		 */
		virtual void multiplyBy( const int& n ) = 0;
		
		/*!	\fn multiplyBy( const AlgebraicInteger* ai  )
		 * 	\brief Performs the multiplication of *this by the parameter
		 * 	\param ai( const AlgebraicInteger* ) Pointer to the multiplicand
		 */
		virtual void multiplyBy( const AlgebraicInteger *ai ) = 0;
		
		/*!	\fn add( const AlgebraicInteger* ai )
		 * 	\brief Performs the addition of *this by the parameter
		 * 	\param ai( const AlgebraicInteger* ) Pointer to the summand
		 */
		virtual void add( const AlgebraicInteger *ai ) = 0;
		
		/*!	\fn substract( const AlgebraicInteger* ai ) = 0
		 * 	\brief Substracts the parameter from *this
		 * 	\param ai( const AlgebraicInteger* ) Pointer to the subtrahend
		 */
		virtual void substract( const AlgebraicInteger *ai ) = 0;
		
		/*!	\fn opp( const AlgebraicInteger* ai )
		 * 	\brief Computes the opp
		 */
		virtual void opp() = 0;

		/*!	\fn bool bIsLessThan( const int& n ) const = 0
		 * 	\brief Check whether *this < n
		 */
		virtual bool bIsLessThan( const int& n ) const = 0;
		
		/*!	\fn bool bIsLessThan( const AlgebraicInteger& ai )const = 0
		 * 	\brief Check whether *this < ai
		 */
		virtual bool bIsLessThan( const AlgebraicInteger& ai ) const = 0;
		
		/*!	\fn bIsLessOEThan( const AlgebraicInteger& ai ) const = 0
		 * 	\brief Check whether *this <= ai
		 */
		virtual bool bIsLessOEThan( const AlgebraicInteger& ai ) const = 0;
		
		/*!	\fn bIsGreaterThan( const int& n ) const = 0
		 * 	\brief Check whether *this > n
		 */
		virtual bool bIsGreaterThan( const int& n ) const = 0;
		
		/*!	\fn bIsGreaterOEThan( const int& n ) const = 0
		 * 	\brief Check whether *this >= n
		 */
		virtual bool bIsGreaterOEThan( const int& n ) const = 0;
		
		/*!	\fn bIsEqualTo( const AlgebraicInteger& n ) const = 0
		 * 	\brief Check whether *this == ai
		 */
		virtual bool bIsEqualTo( const AlgebraicInteger& ai ) const = 0;
		
		/*!	\fn bIsEqualTo( const int& n ) const = 0
		 * 	\brief Check whether *this == n
		 */
		virtual bool bIsEqualTo( const int& n ) const = 0;
		
		/*!	\fn operator==( const AlgebraicInteger& ai ) const
		 * 	\brief Check whether *this == ai
		 */
		bool operator==( const AlgebraicInteger& ai ) const;
		
		/*!	\fn operator!=( const int& n ) const
		 * 	\brief Check whether *this != n
		 */
		bool operator!=( const int& n ) const;
		
		/*!	\fn operator!=( const AlgebraicInteger& ai ) const
		 * 	\brief Check whether *this != ai
		 */
		bool operator!=( const AlgebraicInteger& ai ) const;
		
		friend ostream& operator<<( ostream& , AlgebraicInteger const & );
		
		/*!	\fn ostream& print( ostream& ) const
		 * 	\brief Print the number
		 */
		virtual ostream& print( ostream& ) const;
		
		/*!	\fn string to_string( const string& strFormat = "generic", const bool& bProtect = false ) const = 0
		 * 	\brief Convert to string
		 * 	\param strFormat( const string& ): Format for the output; can be: generic, mathematica, pari, latex
		 * 	\param bProtect( const bool& ): If true, avoid certain characters (for example if the string is used in for a file name)
		 */
		virtual string to_string( const string& strFormat = "generic", const bool& bProtect = false ) const = 0;
		
		/*!	\fn double to_double() const = 0
		 * 	\brief Convert to double (approximation)
		 */
		virtual double to_double() const = 0;
		
		/*!	\fn string get_classname() const = 0
		 * 	\brief Return the type (string)
		 */
		virtual string get_classname() const = 0;
		
	private:
		AlgebraicInteger& operator=( AlgebraicInteger const & ); ///< = is not defined ; use set
};

/*!	\fn isLessThanPtrAlgebraicInteger( AlgebraicInteger* a, AlgebraicInteger* b )
 * 	\brief Check if *a < *b
*/
bool isLessThanPtrAlgebraicInteger( AlgebraicInteger* a, AlgebraicInteger* b );

/*!	\fn isEqualToPtrAlgebraicInteger( AlgebraicInteger* a, AlgebraicInteger* b )
 * 	\brief Check if *a = *b
*/
bool isEqualToPtrAlgebraicInteger( AlgebraicInteger* a, AlgebraicInteger* b );

#endif // ALGEBRAICINTEGER_H
