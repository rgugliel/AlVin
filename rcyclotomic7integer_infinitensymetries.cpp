#include "rcyclotomic7integer_infinitensymetries.h"

RCyclotomic7Integer_InfiniteNSymetries::RCyclotomic7Integer_InfiniteNSymetries( AlVin* alvin )
: InfiniteNSymetries( alvin )
{
	unsigned int i, j;
	
	auto aiVectors( alvin->get_aiVectors() );
	
	riVectorsProducts = vector< vector< RCyclotomic7Integer > >( iVectorsCount, vector< RCyclotomic7Integer >( iVectorsCount, 0 ) );
	iVectorsNormsFloor = vector< long int >( iVectorsCount, 0 );
	
	// -------------------------------------------------
	// Quadratic form
	for( auto ai : aiQF )
		rciQF.push_back( Rational<RCyclotomic7Integer>( *dynamic_cast< RCyclotomic7Integer* >( ai ) ) );
	
	// -------------------------------------------------
	// Vectors
	for( i = 0; i < iVectorsCount; i++ )
	{
		auto v( vector< Rational<RCyclotomic7Integer> >( 0 ) );
		
		Matrix< Rational<RCyclotomic7Integer>, Dynamic, 1 > v2;
		
		for( j = 0; j <= iDimension; j++ )
			v.push_back( Rational<RCyclotomic7Integer>( *dynamic_cast< RCyclotomic7Integer* >( aiVectors[i][j] ) ) );
		
		v2 = Matrix< Rational<RCyclotomic7Integer>, 1, Dynamic >::Map( &v[0], 1, iDimension+1 );
		rciVectorsC.push_back( v2 );
	}
	
	// -------------------------------------------------
	// Products
	for( i = 0; i < iVectorsCount; i++ )
	{	
		AlgebraicInteger* aiNorm( alvin->aiBilinearProduct( aiVectors[i], aiVectors[i] ) );
		riVectorsProducts[i][i] = *dynamic_cast< RCyclotomic7Integer* >( aiNorm );
		
		mpz_class iT( riVectorsProducts[i][i].floor() );
		if( !iT.fits_slong_p() )
			throw( string( "RCyclotomic7Integer_InfiniteNSymetries::RCyclotomic7Integer_InfiniteNSymetries: Number too big" ) );
		
		iVectorsNormsFloor[i] = iT.get_si();
		
		delete aiNorm;
		
		for( j = i+1; j < iVectorsCount; j++ )
		{
			AlgebraicInteger* aiProd( alvin->aiBilinearProduct( aiVectors[i], aiVectors[j] ) );
			riVectorsProducts[i][j] = riVectorsProducts[j][i] =  *dynamic_cast< RCyclotomic7Integer* >( aiProd );
			delete aiProd;
		}
	}
	
	// -------------------------------------------------
	// Weights of dotted lines
	qiDottedWeights = vector< vector< Rational<RCyclotomic7Integer> > >( iVectorsCount, vector< Rational<RCyclotomic7Integer> >( iVectorsCount, 0 ) );
	for( i = 0; i < iVectorsCount; i++ )
	{
		for( j = i + 1; j < iVectorsCount; j++ )
		{
			if( iGraphMatrix[i][j] == 2 ) // dotted
			{
				AlgebraicInteger* aiNum( alvin->aiBilinearProduct( aiVectors[i], aiVectors[j] ) );
				AlgebraicInteger* aiNorm1( alvin->aiBilinearProduct( aiVectors[i], aiVectors[i] ) );
				AlgebraicInteger* aiNorm2( alvin->aiBilinearProduct( aiVectors[j], aiVectors[j] ) );
				
				aiNum->multiplyBy( aiNum );
				aiNorm1->multiplyBy( aiNorm2 );
				
				qiDottedWeights[i][j] = qiDottedWeights[j][i] = Rational<RCyclotomic7Integer>( *dynamic_cast< RCyclotomic7Integer* >( aiNum ), *dynamic_cast< RCyclotomic7Integer* >( aiNorm1 ) );
				
				delete aiNum;
				delete aiNorm1;
				delete aiNorm2;
			}
		}
	}
}

RCyclotomic7Integer_InfiniteNSymetries::~RCyclotomic7Integer_InfiniteNSymetries()
{

}


Rational<RCyclotomic7Integer> RCyclotomic7Integer_InfiniteNSymetries::rciVectorNorm( const Matrix< Rational<RCyclotomic7Integer>, -1, 1 >& rciV ) const
{
	Rational<RCyclotomic7Integer> rriRes( - rciQF[0] * rciV(0,0) * rciV(0,0) );
	
	for( unsigned int j(1); j < iVectorSize; j++ )
		rriRes += rciQF[j] * rciV(j,0)* rciV(j,0);
	
	return rriRes;
}

Rational<RCyclotomic7Integer> RCyclotomic7Integer_InfiniteNSymetries::rciVectorsProduct( const Matrix< Rational<RCyclotomic7Integer>, Dynamic, 1 >& rciV1, const Matrix< Rational<RCyclotomic7Integer>, Dynamic, 1 >& rciV2 ) const
{
	Rational<RCyclotomic7Integer> rciRes( - rciQF[0] * rciV1(0,0) * rciV2(0,0) );
	
	for( unsigned int j(1); j < iVectorSize; j++ )
		rciRes += rciQF[j] * rciV1(j,0)* rciV2(j,0);
	
	return rciRes;
}

vector< GraphInvolution > RCyclotomic7Integer_InfiniteNSymetries::FindIsomorphismsInSubgraph( const vector< unsigned int >& iVertices )
{
	vector< GraphInvolution > iIsomorphisms;
	
	unsigned int iVerticesCount( iVertices.size() ), i, j, k;
	
	// -------------------------------------------------------------
	// Do the corresponding vectors span the whole space?
	Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > mVectorsVertices_temp;
	mVectorsVertices_temp.resize( iVectorSize, iVerticesCount );
	
	for( i = 0; i < iVerticesCount; i++ )
		mVectorsVertices_temp.col(i) = rciVectorsC[iVertices[i]];
	
	FullPivLU<Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic >> lu_mVectorsVertices_temp(mVectorsVertices_temp);
	if( lu_mVectorsVertices_temp.rank() < iVectorSize )
		return iIsomorphisms;
	
	// -------------------------------------------------------------
	// We extract a basis of the space, consisting of a subset of vectors designated by iVertices
	auto iPerm( lu_mVectorsVertices_temp.permutationQ().indices() );

	vector< unsigned int > iLinearlyIndependantVertices( iVectorSize, 0 ); ///< Strictly speaking, the indices i (in iVertices!) such that the set vectors[ iVertices[i] ] is linearly indpendant
	
	std::sort(iPerm.derived().data(), iPerm.derived().data()+iVectorSize);
	Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > mVectorsVerticesBasis;
	mVectorsVerticesBasis.resize( iVectorSize, iVectorSize );
	for( i = 0; i < iVectorSize; i++ )
	{
		iLinearlyIndependantVertices[i] = iPerm[i];
		mVectorsVerticesBasis.col(i) = mVectorsVertices_temp.col(iPerm[i]);
	}
	
	// -------------------------------------------------------------
	// Constructing the graph
	igraph_t graph;
	igraph_matrix_t adjMat;
	igraph_matrix_init( &adjMat, iVerticesCount, iVerticesCount );
	
	// vertices coloring
	igraph_vector_int_t verticesColors;
	igraph_vector_int_init( &verticesColors, iVerticesCount );
	
	for( i = 0; i < iVerticesCount; i++ )
	{
		VECTOR(verticesColors)[i] = iVectorsNormsFloor[iVertices[i]];
		
		for( j = i + 1; j < iVerticesCount; j++ )
			MATRIX( adjMat, i, j ) = iGraphMatrix[iVertices[i]][iVertices[j]];
	}
	
	if( igraph_adjacency( &graph, &adjMat, IGRAPH_ADJ_UNDIRECTED ) != IGRAPH_SUCCESS )
		throw( string( "InfiniteNSymetries::FindInvolutionInSubgraph: Error creating graph" ) );
	
	// -------------------------------------------------------------
	// Finding isomorphisms
	igraph_vector_ptr_t maps;
	igraph_vector_ptr_init( &maps, 1 );
	igraph_get_isomorphisms_vf2( &graph, &graph, &verticesColors, &verticesColors, 0, 0, &maps, 0, 0, 0 );
	long int iIsomorphismsCount( igraph_vector_ptr_size( &maps ) );

	for( int is( 1 ) ; is < iIsomorphismsCount; is++ ) // Don't consider the identity
	{
		bool bIsInvolution( true ), bIsIsomorphism( true );
		vector< unsigned int > iPermutation( iVerticesCount, -1 );
		unsigned int iFixedPoints( 0 );
		
		igraph_vector_t* v( (igraph_vector_t*)igraph_vector_ptr_e( &maps, is ) );
		
		// ---------------------------------------------------
		// Checking that it is an involution and getting the permutation
		for( k = 0; k < iVerticesCount; k++ )
		{
			iPermutation[k] = (unsigned int)VECTOR(*v)[k];
			
			if( iPermutation[k] == k )
				iFixedPoints++;
			
			if( (unsigned int)VECTOR(*v)[(unsigned int)VECTOR(*v)[k]] != k )
			{
				bIsInvolution = false;
				break;
			}
		}
		
		igraph_vector_destroy(v);
		free(v);
		
		if( !bIsInvolution )
			continue;
		
		// ---------------------------------------------------
		// Checking that it is an isomorphism
		for( i = 0; i < iVerticesCount && bIsIsomorphism; i++ )
		{
			for( j = i + 1; j < iVerticesCount; j++ )
			{
				// Weight not preserved or dotted and weight of the dotted not preserved
				if( ( iGraphMatrix[iVertices[i]][iVertices[j]] != iGraphMatrix[iVertices[iPermutation[i]]][iVertices[iPermutation[j]]] )
					|| ( iGraphMatrix[iVertices[i]][iVertices[j]] == 2 && !bDottedSameWeight( iVertices[i], iVertices[j], iVertices[iPermutation[i]], iVertices[iPermutation[j]] ) ) )
				{
					bIsIsomorphism = false;
					break;
				}
				
				// If does not preserve some angle
				if( riVectorsProducts[iVertices[i]][iVertices[j]] != riVectorsProducts[iVertices[iPermutation[i]]][iVertices[iPermutation[j]]] )
				{
					bIsIsomorphism = false;
					break;
				}
			}
		}
		
		if( !bIsIsomorphism ) // Does not preserve the dotted
			continue;
		
		// ---------------------------------------------------
		// Here, we still have to check it preserves the lattice generated by the canonical basis
		Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > iIntegralSymmetry; ///< Permutation matrix wrt to a basis contained in the vectors corresponding to iVertices
		Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > mVectorsVerticesBasisInv( mVectorsVerticesBasis.inverse() );
		
		iIntegralSymmetry.resize( iVectorSize, iVectorSize );
		for( i = 0; i < iVectorSize; i++ )
		{
			iIntegralSymmetry.col(i) = mVectorsVerticesBasisInv * rciVectorsC[iVertices[iPermutation[iLinearlyIndependantVertices[i]]]];
		}
		
		Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > mCanBasis( mVectorsVerticesBasis * iIntegralSymmetry * mVectorsVerticesBasisInv );
		bool bPreservesLattice( true );
		
		for( i = 0; i < iVectorSize && bPreservesLattice; i++ )
		{
			for( j = 0; j < iVectorSize; j++ )
			{
				if( !mCanBasis(i,j).get_hasDenominatorOne() )
				{
					bPreservesLattice = false;
					break;
				}
			}
		}
		
		if( bPreservesLattice )
			iIsomorphisms.push_back( { iVertices, iPermutation, iLinearlyIndependantVertices, (iVerticesCount+iFixedPoints)/2 } );
	}

	// ----------------------------------------
	// Freeing memory
	igraph_destroy( &graph );
	igraph_matrix_destroy( &adjMat );
	
	igraph_vector_destroy( (igraph_vector_t*)igraph_vector_ptr_e( &maps, 0 ) ); // identity
	free( (igraph_vector_t*)igraph_vector_ptr_e( &maps, 0 ) ); // identity
	
	igraph_vector_ptr_destroy( &maps );
	igraph_vector_int_destroy( &verticesColors );
	
	return iIsomorphisms;
}

bool RCyclotomic7Integer_InfiniteNSymetries::FindIntegralSymmetryFromSubgraph( const vector< unsigned int >& iVertices )
{
	unsigned int i, j;
	unsigned int iVerticesCount( iVertices.size() );
	
	// -------------------------------------------------------------
	// Does the corresponding polyhedron have at least one vertex?
	vector< vector< unsigned int > > iCox( vector< vector< unsigned int > >( iVerticesCount, vector< unsigned int >( iVerticesCount, 2 ) ) );
	for( i = 0; i < iVerticesCount; i++ )
	{
		for( j = i + 1; j < iVerticesCount; j++ )
			iCox[i][j] = iCox[j][i] = iCoxeterMatrix[iVertices[i]][iVertices[j]];
	}
	
	CoxIter ci( iCox, iDimension );
	ci.exploreGraph();
	ci.computeGraphsProducts();
	ci.bEulerCharacteristicFVector();
	vector< unsigned int > iFV( ci.get_iFVector() );
	if( !iFV[0] )
		return bFinished;
	
	// -------------------------------------------------------------
	// Finding the isomorphisms
	vector< GraphInvolution > iIsomorphisms( FindIsomorphismsInSubgraph( iVertices ) );
	if( !iIsomorphisms.size() )
		return bFinished;
	
	// -------------------------------------------------------------
	// Going through
	WorkWithIsomorphisms( iVertices, iIsomorphisms );
	
	return bFinished;
}

bool RCyclotomic7Integer_InfiniteNSymetries::bDottedSameWeight( const unsigned int& v1, const unsigned int& v2, const unsigned int& w1, const unsigned int& w2 ) const
{
	return ( qiDottedWeights[v1][v2] == qiDottedWeights[w1][w2] );
}

void RCyclotomic7Integer_InfiniteNSymetries::WorkWithIsomorphisms( const vector< unsigned int >& iVertices, const vector< GraphInvolution >& iIsomorphisms )
{
	unsigned int i, j, k, iVerticesCount( iVertices.size() );
	
	// -------------------------------------------
	// Vectors corresponding to the vertices
	Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > m;
	m.resize( iVectorSize, iVerticesCount );
	
	for( unsigned int i(0); i < iVerticesCount; i++ )
		m.col(i) = rciVectorsC[i];
	
	// --------------------------------------------
	// For each isomorphism
	for( auto iIso : iIsomorphisms )
	{
		// ------------------------------------------------
		// The column of m span the space of fixed points
		Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > m;
		m.resize( iVectorSize, iIso.iOrbitCount );
		
		j = 0;
		for( i = 0; i < iVerticesCount && j < iIso.iOrbitCount; i++ )
		{
			if( iIso.iPermutation[i] == i )
				m.col(j++) = rciVectorsC[iVertices[i]];
			else if( iIso.iPermutation[i] < i )
				m.col(j++) = rciVectorsC[iVertices[i]] + rciVectorsC[ iVertices[iIso.iPermutation[i]]];
		}
		
		// We get a basis of the space of fixed points
		FullPivLU<Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic >> lu(m);
		Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > basisFixedPoints( lu.image(m) );
		
		unsigned int iFixedPointDimTemp( basisFixedPoints.cols() );
		if( iFixedPointDimTemp == iVectorSize ) // Soooo much fixed points
			continue;
		
		#pragma omp critical
		{
			if( !rciBasisFixedPoints.cols() ) // These are the first fixed points
			{
				rciBasisFixedPoints = basisFixedPoints;
				iFixedPointsDimension = iFixedPointDimTemp;
				usefulInvolutions.push_back( iIso );
			}
			else // We already have fixed points
			{
				// The goal here is to take the intersection of the previsous fixed points space and the new one
				Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > tempVectors( rciBasisFixedPoints );
				
				tempVectors.conservativeResize( iVectorSize, tempVectors.cols() + iFixedPointDimTemp );
				for( j = 0; j < iFixedPointDimTemp; j++ )
					tempVectors.col(j + iFixedPointsDimension) = basisFixedPoints.col(j);
				
				FullPivLU<Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > > luTemp(tempVectors);
				Matrix< Rational<RCyclotomic7Integer>, Dynamic, Dynamic > ker( luTemp.kernel() );
				
				if( luTemp.dimensionOfKernel() == 0 )
				{
					bFinished = true;
					iFixedPointsDimension = 0;
					rciBasisFixedPoints.resize( iVectorSize, 0 );
					usefulInvolutions.push_back( iIso );
				}
				
				if( ker.cols() < iFixedPointsDimension && !bFinished ) // If we earned something
				{
					usefulInvolutions.push_back( iIso );
					unsigned int iIntersectionDimension( ker.cols() );
					rciBasisFixedPoints.resize( iVectorSize, iIntersectionDimension );
									
					for( j = 0; j < iIntersectionDimension; j++ )
					{
						Matrix< Rational<RCyclotomic7Integer>, Dynamic, 1 > v( ker(0,j)*tempVectors.col(0) );
						
						for( k = 1; k < iFixedPointsDimension; k++ )
							v += ker(k,j)*tempVectors.col(k);
						
						rciBasisFixedPoints.col(j) = v;
					}
					
					iFixedPointsDimension = ker.cols();
					if( iFixedPointsDimension == 1 )
					{	
						if( rciVectorNorm( rciBasisFixedPoints.col(0) ) >= 0 )
						{
							bFinished = true;
						}
					}
					else if( iFixedPointsDimension == 2 )
					{
						// If the basis has two orthogonal vector of positive norm
						if( rciVectorNorm( rciBasisFixedPoints.col(0) ) >= 0 && rciVectorNorm( rciBasisFixedPoints.col(1) ) >= 0 && rciVectorsProduct( rciBasisFixedPoints.col(0), rciBasisFixedPoints.col(1) ) == 0 )
							bFinished = true;
					}
				}
			}
		}
		
		if( bFinished )
			return;
	}
}

void RCyclotomic7Integer_InfiniteNSymetries::print_basisFixedPoints( const string& strSpacer ) const
{
	for( unsigned int i(0); i < iFixedPointsDimension; i++ )
		cout << strSpacer << Matrix<Rational<RCyclotomic7Integer>, 1, Eigen::Dynamic>( rciBasisFixedPoints.col(i).transpose() ) << endl;
}
