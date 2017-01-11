#include "rcyclotomic7integer_alvin.h"

using namespace RCyclotomic7Integer_constants;

std::string RCyclotomic7Integer_AlVin::get_strField() const
{
	return string( "RC7" );
}

bool RCyclotomic7Integer_AlVin::PreRun()
{
	return true;
}

RCyclotomic7Integer_AlVin::RCyclotomic7Integer_AlVin( const vector< RCyclotomic7Integer >& rciQuadraticFormCoeffs, const string &strOuputMathematicalFormat, const bool& bWriteInfo, const bool& bDebug )
: AlVin( strOuputMathematicalFormat, bWriteInfo, bDebug )
{
	RCyclotomic7Integer::initialize();
	
	bool bNegativeFound( false );
	
	// -------------------------------------------------
	// Checking the coefficients
	for( auto i : rciQuadraticFormCoeffs )
	{
		if( i.bIsLessThan( 0 ) )
		{
			if( !bNegativeFound )
			{
				RCyclotomic7Integer *rci( new RCyclotomic7Integer( { -i.iC[0], -i.iC[1], -i.iC[2] } ) );

				RCyclotomic7Integer rciTemp( *rci );
				
				rciTemp.conjugate( 2 );
				if( rciTemp.bIsGreaterThan( 0 ) )
					throw( string( "Quadratic form is not admissible" ) );
				
				rciTemp.conjugate( 2 );
				if( rciTemp.bIsGreaterThan( 0 ) )
					throw( string( "Quadratic form is not admissible" ) );
				
				aiQF.insert( aiQF.begin(), rci );
				
				bNegativeFound = true;
			}
			else
				throw( string( "Quadratic form has incorrect signature" ) );
		}	
		else if( i.bIsGreaterThan( 0 ) )
		{
			RCyclotomic7Integer *rci( new RCyclotomic7Integer( i ) );
			
			RCyclotomic7Integer rciTemp( *rci );
			
			rciTemp.conjugate( 2 );
			if( rciTemp.bIsLessThan( 0 ) )
				throw( string( "Quadratic form is not admissible" ) );
			
			rciTemp.conjugate( 2 );
			if( rciTemp.bIsLessThan( 0 ) )
				throw( string( "Quadratic form is not admissible" ) );
			
			aiQF.insert( lower_bound( aiQF.begin() + ( bNegativeFound ? 1 : 0 ), aiQF.end(), rci ), rci );
		}
	}
	
	if( !bNegativeFound )
		throw( string( "Quadratic form has incorrect signature" ) );
	
	// -------------------------------------------------
	// global work
	initializations();
	rciVectorCurrent = vector< RCyclotomic7Integer >( iDimension + 1, 0 );
	
	// -------------------------------------------------
	// local copy of the quadratic form (to avoid dynamic_cast)
	RCyclotomic7Integer rci2( 2 ); 
	for( unsigned int i(0); i <= iDimension; i++ )
	{
		RCyclotomic7Integer* rci( dynamic_cast< RCyclotomic7Integer* >( aiQF[i] ) );
		rciQF.push_back( RCyclotomic7Integer( *rci ) );
		
		rci2QF.push_back( RCyclotomic7Integer( *rci ) );
		rci2QF[i].multiplyBy( &rci2 );
	}
	
	// -------------------------------------------------
	// Possible values for (e,e)
	findPossibleNorms2();
	
	// -------------------------------------------------
	// Displaying some information
	if( bWriteInfo )
		print_initialInformation();
}

void RCyclotomic7Integer_AlVin::findVector(AlgebraicInteger* aiX0, AlgebraicInteger* aiNorm2)
{
	RCyclotomic7Integer* rci0( dynamic_cast< RCyclotomic7Integer* >( aiX0 ) );
	RCyclotomic7Integer* rciNorm2( dynamic_cast< RCyclotomic7Integer* >( aiNorm2 ) );
	
	// -----------------------------------------------------
	// Preliminary work
	rciBilinearProducts = vector< RCyclotomic7Integer >( iVectorsCount_second, 0 );
	
	// initial value of iSumComp
	RCyclotomic7Integer rciSumComp( *rciNorm2 );
	RCyclotomic7Integer rciTemp( *rci0 );
	rciTemp.multiplyBy( rci0 );
	rciTemp.multiplyBy( &rciQF[0] );
	rciSumComp.add( &rciTemp );
	
	// initial value of qiBilinearProducts
	for( unsigned int i( 0 ); i < iVectorsCount_second; i++ )
	{
		// = - rci[0] * k0 * rciVectors[i + iDimension][0]
		rciBilinearProducts[i].set( &rciQF[0] );
		rciBilinearProducts[i].multiplyBy( -1 );
		rciBilinearProducts[i].multiplyBy( rci0 );
		rciBilinearProducts[i].multiplyBy( rciVectors[i + iDimension][0] );
	}
	
	rciVectorCurrent[0].set( rci0 );
	
	findVector( rci0, dynamic_cast< RCyclotomic7Integer* >( aiNorm2 ), 1, rciSumComp, *rci0 );
}

void RCyclotomic7Integer_AlVin::print_initialInformationChild() const
{
	
	if( strOuputMathematicalFormat == "mathematica" ) // quadratic form
	{
		cout << "l1 := 2*Cos[2*\\[Pi]/7];\nl2 := 2*Cos[4*\\[Pi]/7];\nl3 := 2*Cos[6*\\[Pi]/7];\nRCI7[l_] := l[[1]] * l1 + l[[2]] * l2 + l[[3]] * l3;" << endl;
	}
	else if( strOuputMathematicalFormat == "pari" )
	{
		cout << "RCI7 = (x,y,z) -> 2*( x*cos(2*Pi/7) + y*cos(4*Pi/7) + z*cos(6*Pi/7) );" << endl;
	}
}


void RCyclotomic7Integer_AlVin::findVector(RCyclotomic7Integer* rci0, RCyclotomic7Integer* rciNorm2, unsigned int iIndex, RCyclotomic7Integer rciSumComp, RCyclotomic7Integer rciGCDComponents)
{
	// ----------------------------------------
	// Frequent values
	static interval gaol_minus1( interval( -1 ) );
	static interval gaol_2( interval( 2 ) );
	static interval gaol_3( interval( 3 ) );
	static interval gaol_7( interval( 7 ) );
	
	if( rciSumComp.iC[0] == 0 && rciSumComp.iC[1] == 0 && rciSumComp.iC[2] == 0 ) // We have a candidate
	{
		for( ; iIndex <= iDimension; iIndex++ )
			rciVectorCurrent[iIndex] = 0;
		
		if( rciGCDComponents.bIsInvertible() )
			addCandidate();
		
		return;
	}
	
	// ----------------------------------------------------
	// Varia
	unsigned int j;
	bool bAdmissible, bReturn( false );
	RCyclotomic7Integer rciTemp, rciPartialNorm, rci, rciLastComponent;
	RCyclotomic7Integer rciSumComp_conj2( rciSumComp ), rciSumComp_conj3( rciSumComp );
	rciSumComp_conj2.conjugate(2);
	rciSumComp_conj3.conjugate(3);
	
	// ----------------------------------------------------
	// Variables for the computations of bounds
	interval gaol_b, gaol_a;
	long int aMin, aMax, bMax, cMin, cMax;
	bool bBoundsUncertain( false );
	
	interval gaol_R1( gaol_SQRTQuotient( rciSumComp, rciQF[iIndex] ) );
	
	rciSumComp.conjugate(2);
	rciQF[iIndex].conjugate(2);
	interval gaol_R2( gaol_SQRTQuotient( rciSumComp, rciQF[iIndex] ) );
	
	rciSumComp.conjugate(2);
	rciQF[iIndex].conjugate(2);
	interval gaol_R3( gaol_SQRTQuotient( rciSumComp, rciQF[iIndex] ) );
	
	rciSumComp.conjugate(2);
	rciQF[iIndex].conjugate(2);
	
	// ----------------------------------------------------
	// Eventual simplifications in nice cases
	bool bR1IsAnInt( gaol_R1.is_an_int() ), bR2IsAnInt( gaol_R2.is_an_int() ), bR3IsAnInt( gaol_R3.is_an_int() );
	bool bRiAreLongInt( bR1IsAnInt && bR2IsAnInt && bR3IsAnInt );
	long int iR1( bR1IsAnInt ? floor( gaol_R1.left() ) : 0 ), iR2( bRiAreLongInt ? floor( gaol_R2.left() ) : 0 ), iR3( bRiAreLongInt ? floor( gaol_R3.left() ) : 0 );
	
	// ----------------------------------------------------
	// Bounds for b
	interval gaol_temp( ( gaol_2 * ( gaol_l1 + gaol_l2 ) + gaol_3 * gaol_l3 ) * gaol_R2 );
	
	// Second term: (l1 + 4 * l2 + 2 * l3 ) ) * R3
	gaol_temp += ( gaol_l1 + gaol_2 * ( gaol_2 * gaol_l2 + gaol_l3 ) ) * gaol_R3;
	
	// Third term (min): R1 * ( 2 * l1 + 3 * l2 + 2 * l3 )
	interval gaol_min( gaol_temp + gaol_R1 * ( gaol_2 * ( gaol_l1 + gaol_l3 ) + gaol_3 * gaol_l2 )  );
	
	// Third term (max)
	interval gaol_max( gaol_temp * gaol_minus1 );
	
	gaol_min = ceil( gaol_min / gaol_7 );
	gaol_max = floor( gaol_max / gaol_7 );
	
	if( !gaol_min.is_an_int() )
		bBoundsUncertain = true;
	
	if( !gaol_max.is_an_int() )
		bBoundsUncertain = true;

	bMax = floor( gaol_max.right() );
	
	for( long int b( floor( gaol_min.left() ) ); b <= bMax && !bReturn; b++ )
	{
		gaol_b = interval( b );
			
		// ------------------------------------
		// Bounds for a
		if( bRiAreLongInt )
		{
			// -----------------------------------------------
			// Min
			long int iC0( 2 * iR1 + 3 * iR2 + 2 * iR3 );
			long int iC1( iR1 + 5 * iR2 + iR3 );
			long int iC2( 7 * b + 4 * iR1 + 6 * iR2 + 4 * iR3 );
			
			if( iC0 == iC1 && iC1 && iC2 && !( iC0 % 7 ) ) // the result is an int
				aMin = -iC0 / 7;
			else
			{
				gaol_min = interval( iC0 ) * gaol_l1 + interval( iC1 ) * gaol_l2 + interval( iC2 ) * gaol_l3;
				gaol_min = ceil( gaol_min / gaol_7 );
				
				if( !gaol_min.is_an_int() )
					throw( string( "RCyclotomic7Integer_AlVin::findVector() : Not an int (bound inf for a)" ) );
				
				aMin = floor( gaol_min.left() );
			}
			
			// -----------------------------------------------
			// Max
			iC0 = 3 * iR2 + 2 * iR3;
			iC1 = 5 * iR2 + iR3;
			iC2 = -7 * b + 6 * iR2 + 4 * iR3;
			
			if( iC0 == iC1 && iC1 && iC2 && !( iC0 % 7 ) ) // the result is an int
				aMax = -iC0 / 7;
			else
			{
				gaol_max = interval( -iC0 ) * gaol_l1 + interval( -iC1 ) * gaol_l2 + interval( -iC2 ) * gaol_l3;
				gaol_max = floor( gaol_max / gaol_7 );
				
				if( !gaol_max.is_an_int() )
					throw( string( "RCyclotomic7Integer_AlVin::findVector() : Not an int (bound max for a)" ) );
				
				aMax = floor( gaol_max.left() );
			}
		}
		else
		{
			gaol_min = gaol_minus1 * ( gaol_2 * ( gaol_l1 + gaol_2 * gaol_l3 ) + gaol_l2 );
			gaol_max = gaol_min;
			
			gaol_min *= ( ( gaol_2 * ( gaol_l1 + gaol_l3 ) + gaol_l2 * gaol_3 ) * gaol_b + gaol_l3 * gaol_R2 - gaol_R3 - gaol_R1 );
			gaol_max *= ( ( gaol_2 * ( gaol_l1 + gaol_l3 ) + gaol_l2 * gaol_3 ) * gaol_b - gaol_l3 * gaol_R2 + gaol_R3 );
			
			gaol_min = ceil( gaol_min / gaol_7 );
			gaol_max = floor( gaol_max / gaol_7 );	
			
			if( !gaol_min.is_an_int() )
				bBoundsUncertain = true;
			
			if( !gaol_max.is_an_int() )
				bBoundsUncertain = true;
			
			aMin = floor( gaol_min.left() );
			aMax = floor( gaol_max.right() );
		}
		
		for( long int a( aMin ); a <= aMax && !bReturn; a++ )
		{
			gaol_a = interval( a );
			
			gaol_max = RCyclotomic7Integer( { 2 * a - 2 * b, -b, a - b } ).to_interval();
			
			if( !bR1IsAnInt )
				gaol_min = gaol_max - gaol_SQRTQuotient( rciSumComp, RCyclotomic7Integer( {-1, -2, -2} ) );
			else
				gaol_min = RCyclotomic7Integer( { iR1 + 2 * a - 2 * b, -b, iR1 + a - b } ).to_interval();
			
			gaol_min = ceil( gaol_min );
			gaol_max = floor( gaol_max );
			
			if( !gaol_min.is_an_int() )
				bBoundsUncertain = true;
			
			if( !gaol_max.is_an_int() )
				bBoundsUncertain = true;

			cMin = floor( gaol_min.left() );
			cMax = floor( gaol_max.right() );
			
			for( long int c( cMax ); c >= cMin; c-- )
			{
				rci.iC = { a, b, c };
				
				// ------------------------------------------------------
				// crystallographic condition
				rciTemp.set( &rci2QF[iIndex] );
				rciTemp.multiplyBy( &rci );
				if( !rciTemp.bIsDivisbleBy( rciNorm2 ) )
					continue;
				
				// ------------------------------------------------------
				// Inequalities 1 of the norm system
				RCyclotomic7Integer rciSquare( rci );
				rciSquare.multiplyBy( &rci );
				rciPartialNorm.set( &rciSquare );
				rciPartialNorm.multiplyBy( &rciQF[iIndex] );
				
				if( bBoundsUncertain ) // In this contexte, we have to check if rci <= R1
				{
					if( !rci.bIsGreaterOEThan( 0 ) )
						continue;
					
					if( !rci.bIsLessOEThan( rciSumComp ) )
						continue;
				}
				
				// ------------------------------------------------------
				// Inequalities 2&3 of the norm system
				rciPartialNorm.conjugate(2);
				if( !rciPartialNorm.bIsLessOEThan( rciSumComp_conj2 ) )
					continue;
				
				rciPartialNorm.conjugate(2);
				if( !rciPartialNorm.bIsLessOEThan( rciSumComp_conj3 ) )
					continue;
				
				rciPartialNorm.conjugate(2);
				
				// ------------------------------------------------------
				// If k_j should be less or equal than k_{j-1}
				if( iComponentLessThan[ iIndex ] && !rci.bIsLessOEThan( rciVectorCurrent[ iComponentLessThan[ iIndex ] ] ) )
					continue;
				
				// ------------------------------------------------------
				// Is the vector admissible?
				bAdmissible = true;

				for( j = 0; j < iVectorsCount_second; j++ )
				{
					// rciTemp = rciQF[iIndex] * ( rci - rciLastComponent ) * rciVectors[j + iDimension][iIndex]
					rciTemp.set( &rci );
					rciTemp.substract( &rciLastComponent );
					rciTemp.multiplyBy( &rciQF[iIndex] );
					rciTemp.multiplyBy( rciVectors[j + iDimension][iIndex] );
					rciBilinearProducts[j].add( &rciTemp );
					
					// TODO: garder en mémoire la valeur max?
					
					if( rciBilinearProducts[j].bIsGreaterThan( 0 ) ) // qiVectorCurrent is not admissible
						bAdmissible = false;
				}
				
				rciLastComponent.set( &rci );
				
				if( !bAdmissible )
					break;
				
				rciVectorCurrent[ iIndex ] = rci;
			
				if( iIndex == iDimension && rciSumComp.bIsEqualTo( rciPartialNorm ) )
				{
					rciGCDComponents.gcd( &rci );
					if( rciGCDComponents.bIsInvertible() )
						addCandidate();
					
					bReturn = true;
					break;
				}
				
				if( iIndex < iDimension )
				{
					RCyclotomic7Integer rciSumSub( rciSumComp );
					rciSumSub.substract( &rciPartialNorm );
					
					rciTemp.set( &rciGCDComponents );
					rciTemp.gcd( &rci );

					findVector( rci0, rciNorm2, iIndex + 1, rciSumSub, rciTemp );
				}
			}
		}
	}
	
	for( unsigned int i( 0 ); i < iVectorsCount_second; i++ )
	{
		rciTemp.set( &rciQF[iIndex] );
		rciTemp.multiplyBy( rciVectors[i + iDimension][iIndex] );
		rciTemp.multiplyBy( &rciLastComponent );
		rciBilinearProducts[i].substract( &rciTemp );
	}
}

void RCyclotomic7Integer_AlVin::addCandidate()
{
	vector< AlgebraicInteger* > aiV;
	
	for( auto i : rciVectorCurrent )
		aiV.push_back( new RCyclotomic7Integer( i ) );
		
	aiVectors_candidates.push_back( aiV );
}

int RCyclotomic7Integer_AlVin::addVector_iFindWeight(AlgebraicInteger* aiNumerator, AlgebraicInteger* aiDenominator)
{
	RCyclotomic7Integer* rciNumerator( dynamic_cast< RCyclotomic7Integer* >( aiNumerator ) );
	RCyclotomic7Integer* rciDenominator( dynamic_cast< RCyclotomic7Integer* >( aiDenominator ) );
	
	if( rciNumerator->iC == array< mpz_class, 3 >( {-1,-2,-2 } ) && rciDenominator->bIsEqualTo( 4 ) )
		return 7;
		
	return -2;
}

void RCyclotomic7Integer_AlVin::findPossibleNorms2()
{
	vector< AlgebraicInteger* > aiPossibleNorms2( { new RCyclotomic7Integer( 1 ) } ); ///< Possible values for (e,e)
	vector< RCyclotomic7Integer > rciFactors; ///< Prime numbers and fundamental unit
	
	// ----------------------------------------------------------
	// Compute for ever coefficient a of the quadratic form the prime numbers which divide a
	for( unsigned int i(0); i <= iDimension; i++ )
	{
		vector< RCyclotomic7Integer > rciTemp( rciQF[i].rciPrimeFactors() );
		rciFactors.insert( rciFactors.end(), rciTemp.begin(), rciTemp.end() );
	}
	
	// Removing duplicates
	sort( rciFactors.begin(), rciFactors.end() );
	rciFactors = vector< RCyclotomic7Integer >( rciFactors.begin(), unique( rciFactors.begin(), rciFactors.end() ) );
	
	// Adding 2
	rciFactors.push_back( 2 );
	
	// Fundamental units
	rciFactors.push_back( RCyclotomic7Integer( {0,0,-1} ) );
	rciFactors.push_back( RCyclotomic7Integer( {0,-1,-1} ) );
	
	unsigned int iFactorsCount( rciFactors.size() );
	if( iFactorsCount > 8 )
		throw( string( "Number of values for (e,e) is too big" ) );
	
	// ----------------------------------------------------------
	// Powers of two
	vector< unsigned int > pow2( { 1 } );
	for( unsigned int i( 1 ); i < iFactorsCount; i++ )
		pow2[i] = pow2[i-1] * 2;
	
	// ----------------------------------------------------------
	// Compute the products
	unsigned int iMax( pow( 2, rciFactors.size() ) - 1 ), j;
	
	for( unsigned int i = 1; i <= iMax; i++ )
	{
		RCyclotomic7Integer* rciProduct( new RCyclotomic7Integer( 1 ) );
		
		for( j = 0; j < iFactorsCount; j++ )
		{
			if( i & pow2[j] )
				rciProduct->multiplyBy( &rciFactors[j] );
		}
		
		// is the product admissible (i.e. positive under galois embeddings)
		rciProduct->conjugate( 2 );
		if( rciProduct->bIsLessThan( 0 ) )
		{
			delete rciProduct;
			continue;
		}
		rciProduct->conjugate( 2 );
		if( rciProduct->bIsLessThan( 0 ) )
		{
			delete rciProduct;
			continue;
		}
		
		rciProduct->conjugate( 2 );
		
		// yes, we keep it
		aiPossibleNorms2.push_back( rciProduct );
	}
	
	// removing duplicates
	sort( aiPossibleNorms2.begin(), aiPossibleNorms2.end(), isLessThanPtrAlgebraicInteger );
	aiPossibleNorms2 = vector< AlgebraicInteger* >( aiPossibleNorms2.begin(), unique( aiPossibleNorms2.begin(), aiPossibleNorms2.end(), isEqualToPtrAlgebraicInteger ) );
	
	vf = new RCyclotomic7Integer_VFs( aiPossibleNorms2, rciQF[0] );
}

void RCyclotomic7Integer_AlVin::addVectorChild(const vector< AlgebraicInteger* >& aiVector)
{
	vector< RCyclotomic7Integer* > rciV;
	for( auto ai : aiVector )
		rciV.push_back( dynamic_cast< RCyclotomic7Integer* >( ai ) );
	
	rciVectors.push_back( rciV );
}
