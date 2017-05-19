#include "alvin.h"

AlVin::AlVin( const string& strOuputMathematicalFormat, const bool& bWriteInfo, const bool& bDebug )
: bComputeInvariantsPolyhedron( false ),
bWriteInfo( bWriteInfo ),
bDebug( bDebug ),
iCreateImage( -1 ),
iVectorsCount_second( 0 ),
vf( nullptr ),
ptrCI( nullptr ),
strOuputMathematicalFormat( strOuputMathematicalFormat )
{
}

AlVin::~AlVin()
{
	for( auto i : aiQF )
		delete i;
	
	for( auto i : iBilinearProducts )
		delete i;
	
	for( auto v : aiVectors )
	{
		for( auto i : v )
			delete i;
	}
	
	if( vf != nullptr )
		delete vf;
	
	for( auto i : aiVectors_candidates )
	{
		for( auto j : i )
			delete j;
	}
	
	if( ptrCI != nullptr )
		delete ptrCI;
}

void AlVin::initializations()
{
	if( aiQF.size() < 3 )
		throw( string( "Quadratic form has not enough coefficients" ) );
	
	ptrCI = nullptr;
	iBilinearProducts.clear();
	iDimension = aiQF.size() - 1;
	
	for( auto i : aiQF )
		i->removeSquareFactors();
	
	sort( aiQF.begin() + 1, aiQF.end(), isLessThanPtrAlgebraicInteger );
	
	// ------------------------------------------------
	// Test: quotients of the positive elements of the QF should not be the square of an invertible element
	for( unsigned int i(1); i <= iDimension; i++ )
	{
		while( i < iDimension && *aiQF[i+1] == *aiQF[i] )
			i++;
		
		for( unsigned int j(i + 1); j <= iDimension; j++ )
		{
			if( aiQF[j]->bIsDivisbleBy( aiQF[i] ) )
			{
				unique_ptr< AlgebraicInteger > aiTest( aiQF[j]->copy() );
				aiTest->divideBy( aiQF[i] );
				
				if( aiTest->bIsSquareOfIvertible() && get_strField() == "RC7" )
				{
					cout << "Check that " << *aiTest << " is not the square of an invertible element" << endl;
					strFinalInformation += "Check that " + aiTest->to_string() + " is not the square of an invertible element";
				}
			}
		}
	}
	
	// ------------------------------------------------
	// Test: GCD
	unique_ptr< AlgebraicInteger > igcd( aiQF[0]->copy() );
	for( unsigned int i( 1 ); i <= iDimension; i++ )
		igcd->gcd( aiQF[i] );
	
	if( !igcd->bIsInvertible() )
		throw( string( "GCD of the coefficients of the quadratic form should be one" ) );
	
	// ------------------------------------------------
	// Other stuff
	strAlgebraicIntegerType = aiQF[0]->get_classname();
	
	findFirstVectors();
}

void AlVin::findFirstVectors()
{
	iVectorsCount = iDimension;
	iVectorsCount_second = 0;
	
	// ------------------------------------------------
	// We count the multiplicity of the coefficients
	unsigned int iCount( 1 );
	
	unique_ptr< AlgebraicInteger > aiLastCoeff( aiQF[1]->copy() );
	
	for( unsigned i( 2 ); i <= iDimension; i++ )
	{
		if( *aiQF[i] != *aiLastCoeff  )
		{
			iQBlocksSize.push_back( iCount );
			iCount = 0;
		}
		
		aiLastCoeff->set( aiQF[i] );
		iCount++;
	}
	
	iQBlocksSize.push_back( iCount );
	
	// ------------------------------------------------
	// Vectors and gram matrix
	iCoxeterMatrix = vector< vector< unsigned int > >( iDimension, vector< unsigned int >( iDimension, 2 ) );
	iComponentLessThan = vector< unsigned int >( iDimension + 1, 0 );
	
	unsigned int iVectorIndex( 1 ), i;
	for( auto iBlockSize: iQBlocksSize )
	{
		for( i = 1; i <= iBlockSize; i++ )
		{
			vector< AlgebraicInteger* > aiV;
			for( unsigned int j( 0 ); j <= iDimension; j++ )
				aiV.push_back( aiQF[0]->aiCopyToInteger( 0 ) );
			
			if( i >= 2 )
				iComponentLessThan[ iVectorIndex ] = iVectorIndex - 1;
			
			if( i == iBlockSize )
				aiV[ iVectorIndex ]->set( -1 );
			else
			{
				aiV[ iVectorIndex ]->set( -1 );
				aiV[ iVectorIndex + 1 ]->set( 1 );

				iCoxeterMatrix[iVectorIndex-1][iVectorIndex] = iCoxeterMatrix[iVectorIndex][iVectorIndex-1] = i < iBlockSize - 1 ? 3 : 4;
			}
			
			aiVectors.push_back( aiV );
			addVectorChild( aiV );
			
			iVectorIndex++;
		}
	}
}

bool AlVin::Run( unsigned int iMinVectors, unsigned int iMaxVectors, bool bLastCheckFV )
{
	if( !PreRun() )
		return false;
	
	unsigned int iCandidatesCount, i, j;

	while( ( !iMaxVectors || iVectorsCount < iMaxVectors ) )
	{
		vector< AlVinFraction* > vfs( vf->getNextAlVinFraction() );
		for( auto frac : vfs )
			findVector( frac->aiX0, frac->aiNorm2 );
		
		// ------------------------------------------------
		// Are the candidates compatible between themselves
		iCandidatesCount = aiVectors_candidates.size();
		
		for( i = 0; i < iCandidatesCount; i++ )
		{
			for( j = i + 1; j < iCandidatesCount; j++ )
			{
				unique_ptr< AlgebraicInteger > aiProd( aiBilinearProduct( aiVectors_candidates[i], aiVectors_candidates[j] ) );
				
				if(  aiProd->bIsGreaterThan( 0 ) ) // This should not happen
					throw( string( "AlVin::Run: Incompatible candidates" ) );
			}
		}
		
		// ------------------------------------------------
		// Adding the vectors
		for( i = 0; i < iCandidatesCount && ( !iMaxVectors || iVectorsCount < iMaxVectors ); i++ )
		{
			if( bWriteInfo )
			{
				printFoundVector( aiVectors_candidates[i], iVectorsCount+1, false );
			}
			
			try
			{
				addVector( aiVectors_candidates[i] );
				addVectorChild( aiVectors_candidates[i] );
			}
			catch( string& strE )
			{
				aiVectors_candidates[i].clear();
				throw( strE );
			}
			
			// So that we don't try do free the memory (will be done via aiVectors)
			aiVectors_candidates[i].clear();

			// --------------------------------------------------
			// Test volume
			if( iMinVectors && iVectorsCount < iMinVectors )
				continue; // No covolume test
			
			ptrCI = new CoxIter( iCoxeterMatrix, iDimension );
			ptrCI->set_bCheckCofiniteness( true );

			if( !ptrCI->bCanBeFiniteCovolume() )
			{
				delete ptrCI;
				ptrCI = nullptr;
				continue;
			}	
			
			if( ( iVectorsCount == iMaxVectors && !bLastCheckFV ) || ptrCI->isFiniteCovolume() == 1 )
			{
				if( bWriteInfo )
				{
					print_finallInformation();
					
					cout << "Algorithm ended\n" << endl;
				}
				
				string strFilename( to_string( iDimension ) + "-" );
				for( j = 0; j <= iDimension; j++ )
					strFilename += aiQF[j]->to_string( "filename" ) + ( j < iDimension ? "," : "" );

				if( !ptrCI->bWriteGraph( "output/" + strFilename ) || !ptrCI->bWriteGraphToDraw( "output/" + strFilename ) )
					cout << "Error:\n\tCheck that the folder 'output/' exists and is writable" << endl;
				else if( bWriteInfo )
				{
					cout << "Graph written in the file: \n\t" << "output/" + strFilename << ".coxiter" << endl;
					
					if( bWriteInfo )
					{
						if( iCreateImage == 1 || ( iCreateImage == -1 && iVectorsCount <= 25 ) )
						{
							FILE *fin;
							string strCmd( "dot -Tjpg -ooutput/" + strFilename + ".jpg  output/" + strFilename + ".graphviz" );
							if( ( fin = popen( strCmd.c_str(), "r" ) ) )
								pclose( fin );
						}
						else
						{
							cout << "\nCommand to create the image: \n\tdot -Tjpg -ooutput/" << strFilename << ".jpg  output/" << strFilename << ".graphviz\n" << endl;
						}
					}
				}
				
				return true;
			}
			
			if( ptrCI != nullptr )
			{
				delete ptrCI;
				ptrCI = nullptr;
			}
		}
		
		aiVectors_candidates.clear();
	}
	
	CoxIter ci( iCoxeterMatrix, iDimension );
	
	string strFilename( "nt-" + to_string( iDimension ) + "-" );
	for( j = 0; j <= iDimension; j++ )
		strFilename += aiQF[j]->to_string() + ( j < iDimension ? "," : "" );

	ci.bWriteGraph( "output/" + strFilename );
	ci.bWriteGraphToDraw( "output/" + strFilename );
	
	if( bWriteInfo )
	{
		if( iCreateImage == 1 || ( iCreateImage == -1 && iVectorsCount <= 25 ) )
		{
			FILE *fin;
			string strCmd( "dot -Tjpg -ooutput/" + strFilename + ".jpg output/" + strFilename + ".graphviz" );
			if( ( fin = popen( strCmd.c_str(), "r" ) ) )
				pclose( fin );
		}
		else
			cout << "\nCommand to create the image: \n\tdot -Tjpg -ooutput/" << strFilename << ".jpg output/" << strFilename << ".graphviz\n" << endl;
		
		print_finallInformation();
	}
	
	return false;
}

void AlVin::addVector( const vector< AlgebraicInteger* >& iVect )
{
	int iWeight;
	
	AlgebraicInteger* aiNorm( aiBilinearProduct( iVect, iVect ) );
	
	iCoxeterMatrix.push_back( vector< unsigned int >( iVectorsCount + 1, 2 ) );
	
	aiVectors.push_back( iVect );
	for( unsigned int i( 0 ); i < iVectorsCount; i++ )
	{
		AlgebraicInteger* aiProd( aiBilinearProduct( aiVectors[i], iVect ) );

		iWeight = -2;
		
		if( aiProd->bIsEqualTo( 0 ) )
			iWeight = 2;
		else
		{
			AlgebraicInteger* aiNumerator( aiProd->copy() );
			aiNumerator->multiplyBy( aiProd );
			
			AlgebraicInteger* aiDenominator( aiBilinearProduct( aiVectors[i], aiVectors[i] ) );
			aiDenominator->multiplyBy( aiNorm );
			
			if( aiDenominator->bIsLessThan( *aiNumerator ) )
				iWeight = 1; // Dotted edge
			else
			{
				AlgebraicInteger* aiGCD( aiNumerator->copy() );
				aiGCD->gcd( aiDenominator );
				
				aiNumerator->divideBy( aiGCD );
				aiDenominator->divideBy( aiGCD );
				delete aiGCD;
				
				if( aiNumerator->bIsEqualTo( 1 ) )
				{
					if( aiDenominator->bIsEqualTo( 4 ) )
						iWeight = 3;
					else if( aiDenominator->bIsEqualTo( 2 ) )
						iWeight = 4;
					else if( aiDenominator->bIsEqualTo( 1 ) )
						iWeight = 0; // Infty
				}
				else if( aiNumerator->bIsEqualTo( 3 ) )
				{
					if( aiDenominator->bIsEqualTo( 4 ) )
						iWeight = 6;
				}
				
				if( iWeight == -2 )
					iWeight = addVector_iFindWeight( aiNumerator, aiDenominator );
			}
			
			if( iWeight == -2 )
				throw( string( "Unknown weight betweens vectors " + to_string( i + 1 ) + " and " + to_string( iVectorsCount + 1 ) ) );
			
			delete aiDenominator;
			delete aiNumerator;
		}
		
		delete aiProd;
		
		iCoxeterMatrix[i].push_back( iWeight );
		iCoxeterMatrix[iVectorsCount][i] = iWeight;
	}
	
	iVectorsCount++;
	iVectorsCount_second++;
	
	delete aiNorm;
}

void AlVin::printFoundVector( std::vector< AlgebraicInteger* > aiV, const unsigned int& iIndex, const bool& bFirst ) const
{
	if( strOuputMathematicalFormat == "latex" )
	{
		cout << "\te_" << ( iIndex > 9 ? "{" : "" ) << iIndex << ( iIndex > 9 ? "}" : "" ) << " = \\left(";
		for( auto it( aiV.begin() ); it != aiV.end(); ++it )
			cout << ( it == aiV.begin() ? "" : ", " ) << (*it)->to_string( "latex" );
		cout << "\\right)" << endl;
	}
	else if( strOuputMathematicalFormat == "mathematica" )
	{
		if( bFirst )
		{
			cout << "vectors := {{";
			for( auto it( aiV.begin() ); it != aiV.end(); ++it )
				cout << ( it == aiV.begin() ? "" : ", " ) << (*it)->to_string( "mathematica" );
			cout << "}};" << endl;
		}
		else
		{
			cout << "AppendTo[vectors, {";
			for( auto it( aiV.begin() ); it != aiV.end(); ++it )
				cout << ( it == aiV.begin() ? "" : ", " ) << (*it)->to_string( "mathematica" );
			cout << "}]; (* vector " << iIndex << " *); " << endl;
		}
	}
	else if( strOuputMathematicalFormat == "pari" ) 
	{
		if( bFirst )
			cout << "vectors = Mat( [";
		else
			cout << "vectors = concat(vectors, [";
		
		for( auto it( aiV.begin() ); it != aiV.end(); ++it )
			cout << ( it == aiV.begin() ? "" : ", " ) << (*it)->to_string( "pari" );
		cout << "] );" << endl;
	}
	else
	{
		cout << "\te" << iIndex << " = (";
		for( auto it( aiV.begin() ); it != aiV.end(); ++it )
			cout << ( it == aiV.begin() ? "" : ", " ) << **it;
		cout << ")" << endl;
	}
}

void AlVin::print_iQF() const
{
	cout << "Quadratic form (" << iDimension << ",1): ";
	for( auto it( aiQF.begin() ); it != aiQF.end(); it++ )
	{
		if( it == aiQF.begin() )
		{
			unique_ptr< AlgebraicInteger > ptr( (*it)->copy() ); 
			ptr->opp();
			cout << *ptr;
		}
		else
			cout << ", " << **it;
	}
	cout << endl;
	
	cout << "Field of definition: " << get_strField() << endl;
}

void AlVin::print_initialInformationChild() const
{

}

void AlVin::print_initialInformation() const
{
	print_iQF();
	
	if( bDebug )
	{
		cout << "Conditions on the coefficients: \n";
		for( unsigned int i( 1 ); i <= iDimension; i++ )
		{
			if( iComponentLessThan[i] )
				cout << "\tx" << i << " <= x" << iComponentLessThan[i] << endl;
		}
		cout << endl;
		
		const vector< AlgebraicInteger*>* ptraiPN2( vf->get_ptraiPossibleNorm2() );
		unsigned int iMax( ptraiPN2->size() );
		
		cout << "Possible values for (e,e): ";
		for( unsigned int i(0); i < iMax; i++ )
			cout << ( i ? ", " : "" ) << *(*ptraiPN2)[i];
		cout << endl;
	}
	
	cout << endl;
	
	print_initialInformationChild();
	if( strOuputMathematicalFormat == "mathematica" ) // quadratic form
	{
		cout << "F[x_, y_] := Sum[x[[i]]*y[[i]]*qform[[i]], {i, 1, Length[x]}];\nS[x_, y_] := F[x, y]/Sqrt[F[x, x]*F[y, y]];" << endl;
		cout << "qform := {";
		for( unsigned int i(0); i <= iDimension; i++ )
			cout << ( i ? ", " : "-(" ) << aiQF[i]->to_string( "mathematica" ) << ( i ? "" : ")" );
		cout << "};" << endl;
	}
	else if( strOuputMathematicalFormat == "pari" )
	{
		cout << "qform = matdiagonal([";
		for( unsigned int i(0); i <= iDimension; i++ )
			cout << ( i ? ", " : "-(" ) << aiQF[i]->to_string( "pari" ) << ( i ? "" : ")" );
		cout << "]);" << endl;
		
		cout << "S = (x,y) -> qfeval(qform,x,y)/sqrt( qfeval(qform,x) * qfeval(qform,y) );" << endl;
	}
	else
		cout << "Vectors: " << endl;
	
	for( unsigned int i(0); i < iVectorsCount; i++ )
		printFoundVector( aiVectors[i], i+1, i == 0 ); 
}

void AlVin::print_finallInformation() const
{
	if( strOuputMathematicalFormat == "mathematica" ) // quadratic form
	{
		cout << "grammatrix := Table[S[vectors[[i]], vectors[[j]]], {i, 1, Length[vectors]}, {j, 1, Length[vectors]}];" << endl;
	}
	else if( strOuputMathematicalFormat == "pari" )
	{
		cout << "grammatrix = matrix( " << iVectorsCount << "," << iVectorsCount << ", i, j, S(vectors[i,], vectors[j,]) );" << endl;
	}
}

void AlVin::print_vectors() const
{
	cout << "Vectors found: " << endl;
	
	for( unsigned int i( 0 ); i < iVectorsCount; i++ )
	{
		cout << "\te_" << ( i + 1 ) << " = {";
		for( unsigned int j( 0 ); j <= iDimension; j++ )
			cout << ( j ? "," : "" ) << *aiVectors[i][j];
		cout << "}" << endl;
	}
}

vector< vector< unsigned int > > AlVin::get_iCoxeterMatrix() const
{
	return iCoxeterMatrix;
}

unsigned int AlVin::get_iDimension() const
{
	return iDimension;
}

CoxIter* AlVin::get_ptrCI() const
{
	return ptrCI;
}

string AlVin::get_strCoxeterMatrix() const
{
	string strCox;
	
	for( vector< vector< unsigned int > >::const_iterator itRow( iCoxeterMatrix.begin() ); itRow != iCoxeterMatrix.end(); ++itRow )
	{
		if( itRow != iCoxeterMatrix.begin() )
			strCox += "\n";
		
		for( vector< unsigned int >::const_iterator itCol( itRow->begin() ); itCol != itRow->end(); ++itCol )
			strCox += ( itCol == itRow->begin() ? "" : "," ) + to_string( *itCol );
	}
	
	return strCox;
}

string AlVin::get_strQF( const string& strSeparator ) const
{
	string strQF;
	
	for( auto it( aiQF.begin() ); it != aiQF.end(); it++ )
	{
		if( it == aiQF.begin() )
		{
			unique_ptr< AlgebraicInteger > ptr( (*it)->copy() ); 
			ptr->opp();
			strQF = ptr->to_string();
		}
		else
			strQF += strSeparator + (*it)->to_string();
	}
	
	return strQF;
}

string AlVin::get_strAlgebraicIntegerType() const
{
	return strAlgebraicIntegerType;
}

vector< std::vector< AlgebraicInteger* > > AlVin::get_aiVectors() const
{
	return aiVectors;
}

unsigned int AlVin::get_iVectorsCount() const
{
	return iVectorsCount;
}

vector< AlgebraicInteger* > AlVin::get_aiQF() const
{
	return aiQF;
}

std::vector< AlgebraicInteger* > AlVin::get_aiPossibleNorm2() const
{
	return vf->get_aiPossibleNorm2();
}

const std::vector< AlgebraicInteger* >* AlVin::get_ptraiPossibleNorm2() const
{
	return vf->get_ptraiPossibleNorm2();
}

void AlVin::set_bComputeInvariantsPolyhedron(const bool& bValue)
{
	bComputeInvariantsPolyhedron = bValue;
}

void AlVin::set_iCreateImage(const int& iValue)
{
	iCreateImage = iValue;
}

string AlVin::get_strFinalInformation() const
{
	return strFinalInformation;
}
