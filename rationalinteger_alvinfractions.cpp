#include "rationalinteger_alvinfractions.h"

RationalInteger_VFs::RationalInteger_VFs( vector< AlgebraicInteger* > aiPossibleNorms2 )
: AlVinFractions( aiPossibleNorms2 ), iPossibleNorms2_max( dynamic_cast< RationalInteger* >( aiPossibleNorms2[ aiPossibleNorms2.size() - 1 ] )->get_iValue() )
{
}

RationalInteger_VFs::~RationalInteger_VFs()
{

}

void RationalInteger_VFs::computeNextAlVinFractions()
{
	unsigned int iMin( iLastMaximum );
	unsigned int iNumerator;
	
	iLastMaximum = ( iMin + iBatchSize ) * ( iMin + iBatchSize );
	
	// We generate fractions f such that: iMin < f <= iLastMaximum
	
	for( auto aiE : aiPossibleNorms2 )
	{
		unsigned int iE( dynamic_cast< RationalInteger* >( aiE )->get_iValue() );
		
		// ceil( isqrt( x ) ) = isqrt( x - 1 ) + 1
		unsigned int iXMin( iMin ? iSQRT( iE * iMin - 1 ) + 1 : 1 ), iXMax( iSQRT( iE * iLastMaximum ) );
		
		if( iXMin * iXMin == iE * iMin ) // because we want iMin < f and not iMin <= f 
			iXMin++;
		
		
		for( ; iXMin <= iXMax; iXMin++ )
		{
			iNumerator = iXMin * iXMin * iPossibleNorms2_max / iE;
			
			AlVinFraction* vf = new AlVinFraction( new RationalInteger( iXMin ), new RationalInteger( iE ), new RationalInteger( iNumerator ) );
			
			alvinfractions.insert( lower_bound( alvinfractions.begin(), alvinfractions.end(), vf, isLessThanPtrAlVinFraction ), vf );
		}
	}
	
	alvinfractions_it = alvinfractions.begin();
}
