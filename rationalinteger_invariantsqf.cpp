#include "rationalinteger_invariantsqf.h"

InvariantsQF::InvariantsQF( const vector< int >& iQuadraticFormCoeffs )
{
	bool bNegativeFound( false );
	
	for( auto i : iQuadraticFormCoeffs )
	{
		if( i < 0 && !bNegativeFound )
		{
			iQF.insert( iQF.begin(), iRemoveSquareFactors( -i ) );
			bNegativeFound = true;
		}
		else if( i < 0 )
			throw( string( "QUADRATIC_FORM_WRONG_SIGNATURE" ) );
		else if( i > 0 )
			iQF.push_back( iRemoveSquareFactors( i ) );
	}
	
	if( !bNegativeFound )
		throw( string( "QUADRATIC_FORM_WRONG_SIGNATURE" ) );	
	
	sort( iQF.begin() + 1, iQF.end() );
	iDimension = iQF.size() - 1;
	
	// -------------------------
	// Let's compute
	
	computeInvariants();
}

void InvariantsQF::computeInvariants()
{
	// -------------------------------------------------------------------------------------------
	// Decomposition of the coefficients of the quadratic form in product of prime numbers
	primeNumbers();
	
	// -------------------------------------------------------------------------------------------
	// Ramification of B
	computeRamification();
	
	// -------------------------------------------------------------------------------------------
	// Final steps
	strInvariant = "{Q,";
	
	// -------------------------------------------------------------------------------------------
	// iDimension is odd --> more computations
	if( iDimension % 2 && iSignedDeterminant != 1 )
	{
		// We look for the primes which split in Q[ sqrt[ iSignedDeterminant ] ] and are in the ramification set
		vector< string > strSplitPrimes;
		
		for( auto iP : iRamification )
		{
			if( !iP || !( iSignedDeterminant % iP ) )
				continue;
			
			if( iP == 2 )
			{
				if( iSignedDeterminant % 8 == 1 || iSignedDeterminant % 8 == -7 )
					strSplitPrimes.push_back( to_string( 2 ) );
			}
			else
			{
				if( iSignedDeterminant % iP && iJacobiSymbol( iSignedDeterminant, iP ) == 1 )
					strSplitPrimes.push_back( to_string( iP ) );
			}
		}
		
		strInvariant += to_string( iSignedDeterminant ) + ",";
		strInvariant += "{" + implode( ",", strSplitPrimes ) + "}";
	}
	else
	{
		string strRam( implode( ",", iRamification ) );
		str_replace( strRam, "0", "infinity" );
		
		strInvariant += "{" + strRam + "}";
	}
	
	strInvariant += "}";
}

void InvariantsQF::computeRamification()
{
	bool bSign( false );
	
	// -------------------------------------------------------------------------------------------
	// Ramification of s(V)
	if( iQF[0] == 1 )
	{
		iQF_primes[0].push_back( 1 ); // 1 is not really a prime number but we need something
	}

	for( unsigned int i( 0 ); i < iDimension; i++ )
	{
		for( unsigned int j( i + 1 ); j <= iDimension; j++ )
		{
			bSign = i != 0;
			
			for( auto iP1 : iQF_primes[i] )
			{
				for( auto iP2 : iQF_primes[j] )
				{
					if( !bSign )
					{
						if( iP1 != iP2 ) // Ram ( -p, p ) == emptyset
							iRamification = iRamificationProduct( iRamification, iRamificationElementaryQuaternionAlgebra( -iP1, iP2 ) );
					}
					else
					{
						if( iP1 == iP2 ) // (p,p) = (-1,p)
							iRamification = iRamificationProduct( iRamification, iRamificationElementaryQuaternionAlgebra( -1, iP2 ) );
						else
							iRamification = iRamificationProduct( iRamification, iRamificationProduct( iRamificationElementaryQuaternionAlgebra( -iP1, iP2 ), iRamificationElementaryQuaternionAlgebra( -1, iP2 ) ) );
					}
				}
				
				bSign = true;
			}
		}
	}
	
	// -------------------------------------------------------------------------------------------
	// Determinant and discriminant
	iDeterminant = 1;
	for( auto i : iQF )
		iDeterminant = iRemoveSquareFactors( iDeterminant * i );
	
	iPrimesDeterminant = iPrimeFactors< unsigned >( iDeterminant );
	iDeterminant *= -1;
	
	iSignedDeterminant = ( iDimension * ( iDimension + 1 ) ) / 2 % 2 ? -iDeterminant : iDeterminant;
	
	// -------------------------------------------------------------------------------------------
	// Coefficient between s(V) and c(V)
	unsigned int iDimensionM8( iDimension % 8 );
	if( iDimensionM8 == 2 || iDimensionM8 == 3 ) // c(V) = s(V) * ( -1, -d(V) )
	{
		for( auto iP : iPrimesDeterminant )
			iRamification = iRamificationProduct( iRamification, iRamificationElementaryQuaternionAlgebra( -1, iP ) );
	}
	else if( iDimensionM8 == 4 || iDimensionM8 == 5 ) // c(V) = s(V) * ( -1, -1 )
	{
		iRamification = iRamificationProduct( iRamification, vector< unsigned int >{0,2} );
	}
	else if( iDimensionM8 == 6 || iDimensionM8 == 7 ) // // c(V) = s(V) * ( -1, d(V) )
	{
		bSign = false;
		if( !iPrimesDeterminant.size() )
			iPrimesDeterminant.push_back( 1 );

		for( auto iP : iPrimesDeterminant )
		{
			if( bSign )
				iRamification = iRamificationProduct( iRamification, iRamificationElementaryQuaternionAlgebra( -1, iP ) );
			else
			{
				iRamification = iRamificationProduct( iRamification, iRamificationElementaryQuaternionAlgebra( -1, -iP ) );
				bSign = true;
			}
		}
	}
}

vector< unsigned int > InvariantsQF::iRamificationProduct( const vector< unsigned int >& iR1, const vector< unsigned int >& iR2 )
{
	if( iR1 == iR2 )
		return vector< unsigned int >( 0 );
	
	if( !iR1.size() )
		return iR2;
	if( !iR2.size() )
		return iR1;
	
	vector< unsigned int > iTemp, iIntersection, iRamification;

	set_union( iR1.begin(), iR1.end(), iR2.begin(), iR2.end(), inserter(iTemp, iTemp.begin()) );
	set_intersection( iR1.begin(), iR1.end(), iR2.begin(), iR2.end(), inserter( iIntersection, iIntersection.begin() ) );
	
	iTemp = vector< unsigned int >( iTemp.begin(), unique( iTemp.begin(), iTemp.end() ) );
	set_difference( iTemp.begin(), iTemp.end(), iIntersection.begin(), iIntersection.end(), inserter( iRamification, iRamification.begin() ) );
	
	return iRamification;
}

vector< unsigned int > InvariantsQF::iRamificationElementaryQuaternionAlgebra( int a, int b ) const
{
	vector< unsigned int > iRam;
	
	int iA( abs( a ) < abs( b ) ? a : b ), iB( abs( a ) < abs( b ) ? b : a );
	unsigned int iB8( abs( iB ) % 8 ), iAP( iA < 0 ? -iA : iA ), iBP( iB < 0 ? -iB : iB );
	
	if( a == -b )
		return iRam;
	else if( ( iA == -1 && iB == 2 ) )
		return iRam;
	else if( ( iA == -2 && iB == -1 ) || ( iA == -1 && iB == -1 ) )
		return vector< unsigned int >{ 0, 2 };
	else if( iA == -1 && iB > 0 )
		return ( iB8 == 1 || iB8 == 5 ? iRam : vector< unsigned int >{ 2, iBP } );
	else if( iA == -1 && iB < 0 )
		return ( iB8 == 1 || iB8 == 5 ? vector< unsigned int >{ 0, 2 } : vector< unsigned int >{ 0, iBP } );
	else if( iA == 2 && iB < 0 )
		return ( iB8 == 1 || iB8 == 7 ? iRam : vector< unsigned int >{ 2, iBP } );
	else if( iA == -2 && iB > 0 )
		return ( iB8 == 1 || iB8 == 3 ? iRam : vector< unsigned int >{ 2, iBP } );
	else if( iA < 0 && iB > 0 ) // here, two different odd primes
	{
		int iLegendre( iJacobiSymbol( iAP, iBP ) );
		
		if( iBP % 4 == 1 )
			return ( iLegendre == 1 ? iRam : vector< unsigned int >{ iAP, iBP } );
		else // iBP % 4 == 3
		{
			if( iAP % 4 == 1 )
				return ( iLegendre == 1 ? vector< unsigned int >{ 2, iBP } : vector< unsigned int >{ 2, iAP } );
			else // iAP % 4 == 3
				return ( iLegendre == 1 ? vector< unsigned int >{ iAP, iBP } : iRam );
		}
	}
	else if( iB < 0 && iA > 0 ) // here, two different odd primes
	{
		int iLegendre( iJacobiSymbol( iBP, iAP ) );
		
		if( iAP % 4 == 1 )
			return ( iLegendre == 1 ? iRam : vector< unsigned int >{ iAP, iBP } );
		else // iAP % 4 == 3
		{
			if( iBP % 4 == 1 ) 
				return ( iLegendre == 1 ? vector< unsigned int >{ 2, iAP } : vector< unsigned int >{ 2, iBP } );
			else
				return ( iLegendre == 1 ? vector< unsigned int >{ iAP, iBP } : iRam );
		}
	}
	
	throw( string( "UNSPECIFIED CASE: (" + to_string( a ) + "," + to_string( b ) + ") / (" + to_string( iA ) + "," + to_string( iB ) + ")" ) );
	
	return iRam;
}

void InvariantsQF::primeNumbers()
{
	iQF_primes = vector< vector< unsigned int > >( iDimension + 1, vector< unsigned int >( 0 ) );
	
	for( unsigned int i( 0 ); i <= iDimension; i++ )
	{
		if( iQF[i] == 1 )
			continue;
		
		if( i && iQF[i] == iQF[i-1] )
		{
			iQF_primes[i] = iQF_primes[i-1];
			continue;
		}
		
		iQF_primes[i] = iPrimeFactorsWithoutSquares( iQF[i] );
	}
}

vector< unsigned int > InvariantsQF::get_iRamification() const
{
	return iRamification;
}

string InvariantsQF::get_strInvariant() const
{
	return strInvariant;
}
