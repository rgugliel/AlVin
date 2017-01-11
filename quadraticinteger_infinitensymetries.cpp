#include "quadraticinteger_infinitensymetries.h"

QuadraticInteger_InfiniteNSymetries::QuadraticInteger_InfiniteNSymetries( AlVin* alvin )
: InfiniteNSymetries( alvin )
{
	unsigned int i, j;
	
	QuadraticIntegerBig::set_d( QuadraticInteger::d );
	
	auto aiVectors( alvin->get_aiVectors() );
	
	riVectorsProducts = vector< vector< QuadraticIntegerBig > >( iVectorsCount, vector< QuadraticIntegerBig >( iVectorsCount, 0 ) );
	iVectorsNormsFloor = vector< long int >( iVectorsCount, 0 );
	
	// -------------------------------------------------
	// Quadratic form
	for( auto ai : aiQF )
		rqiQF.push_back( Rational<QuadraticIntegerBig>( *dynamic_cast< QuadraticInteger* >( ai ) ) );
	
	// -------------------------------------------------
	// Vectors
	for( i = 0; i < iVectorsCount; i++ )
	{
		auto v( vector< Rational<QuadraticIntegerBig> >( 0 ) );
		
		Matrix< Rational<QuadraticIntegerBig>, Dynamic, 1 > v2;
		
		for( j = 0; j <= iDimension; j++ )
			v.push_back( Rational<QuadraticIntegerBig>( *dynamic_cast< QuadraticInteger* >( aiVectors[i][j] ) ) );
		
		v2 = Matrix< Rational<QuadraticIntegerBig>, 1, Dynamic >::Map( &v[0], 1, iDimension+1 );
		rqiVectorsC.push_back( v2 );
	}
	
	// -------------------------------------------------
	// Products
	for( i = 0; i < iVectorsCount; i++ )
	{	
		AlgebraicInteger* aiNorm( alvin->aiBilinearProduct( aiVectors[i], aiVectors[i] ) );
		riVectorsProducts[i][i] = *dynamic_cast< QuadraticInteger* >( aiNorm );
		
		mpz_class iT( riVectorsProducts[i][i].floor() );
		if( !iT.fits_slong_p() )
			throw( string( "QuadraticInteger_InfiniteNSymetries::QuadraticInteger_InfiniteNSymetries: Number too big" ) );
		
		iVectorsNormsFloor[i] = iT.get_si();
		
		delete aiNorm;
		
		for( j = i+1; j < iVectorsCount; j++ )
		{
			AlgebraicInteger* aiProd( alvin->aiBilinearProduct( aiVectors[i], aiVectors[j] ) );
			riVectorsProducts[i][j] = riVectorsProducts[j][i] =  *dynamic_cast< QuadraticInteger* >( aiProd );
			delete aiProd;
		}
	}
	
	// -------------------------------------------------
	// Weights of dotted lines
	qiDottedWeights = vector< vector< Rational<QuadraticIntegerBig> > >( iVectorsCount, vector< Rational<QuadraticIntegerBig> >( iVectorsCount, 0 ) );
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
				
				qiDottedWeights[i][j] = qiDottedWeights[j][i] = Rational<QuadraticIntegerBig>( *dynamic_cast< QuadraticInteger* >( aiNum ), *dynamic_cast< QuadraticInteger* >( aiNorm1 ) );
				
				delete aiNum;
				delete aiNorm1;
				delete aiNorm2;
			}
		}
	}
}

QuadraticInteger_InfiniteNSymetries::~QuadraticInteger_InfiniteNSymetries()
{

}


Rational<QuadraticIntegerBig> QuadraticInteger_InfiniteNSymetries::rqiVectorNorm( const Matrix< Rational<QuadraticIntegerBig>, -1, 1 >& rqiV ) const
{
	if( rqiV.rows() != iVectorSize )
		throw( string( "QuadraticInteger_InfiniteNSymetries::rqiVectorNorm: Invalid size" ) );
	
	Rational<QuadraticIntegerBig> rriRes( - rqiQF[0] * rqiV(0,0) * rqiV(0,0) );
	
	for( unsigned int j(1); j < iVectorSize; j++ )
		rriRes += rqiQF[j] * rqiV(j,0)* rqiV(j,0);
	
	return rriRes;
}

Rational< QuadraticIntegerBig > QuadraticInteger_InfiniteNSymetries::rqiVectorsProduct( const Matrix< Rational<QuadraticIntegerBig>, Dynamic, 1 >& rqiV1, const Matrix< Rational<QuadraticIntegerBig>, Dynamic, 1 >& rqiV2 ) const
{
	Rational<QuadraticIntegerBig> rqiRes( - rqiQF[0] * rqiV1(0,0) * rqiV2(0,0) );
	
	for( unsigned int j(1); j < iVectorSize; j++ )
		rqiRes += rqiQF[j] * rqiV1(j,0)* rqiV2(j,0);
	
	return rqiRes;
}


vector< GraphInvolution > QuadraticInteger_InfiniteNSymetries::FindIsomorphismsInSubgraph( const vector< unsigned int >& iVertices )
{
	vector< GraphInvolution > iIsomorphisms;
	
	unsigned int iVerticesCount( iVertices.size() ), i, j, k;
	
	// -------------------------------------------------------------
	// Do the corresponding vectors span the whole space?
	Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > mVectorsVertices_temp;
	mVectorsVertices_temp.resize( iVectorSize, iVerticesCount );
	
	for( i = 0; i < iVerticesCount; i++ )
		mVectorsVertices_temp.col(i) = rqiVectorsC[iVertices[i]];
	
	FullPivLU<Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic >> lu_mVectorsVertices_temp(mVectorsVertices_temp);
	if( lu_mVectorsVertices_temp.rank() < iVectorSize )
		return iIsomorphisms;
	
	// -------------------------------------------------------------
	// We extract a basis of the space, consisting of a subset of vectors designated by iVertices
	auto iPerm( lu_mVectorsVertices_temp.permutationQ().indices() );

	vector< unsigned int > iLinearlyIndependantVertices( iVectorSize, 0 ); ///< Strictly speaking, the indices i (in iVertices!) such that the set vectors[ iVertices[i] ] is linearly indpendant
	
	std::sort(iPerm.derived().data(), iPerm.derived().data()+iVectorSize);
	Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > mVectorsVerticesBasis;
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
		Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > iIntegralSymmetry; ///< Permutation matrix wrt to a basis contained in the vectors corresponding to iVertices
		Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > mVectorsVerticesBasisInv( mVectorsVerticesBasis.inverse() );
		
		iIntegralSymmetry.resize( iVectorSize, iVectorSize );
		for( i = 0; i < iVectorSize; i++ )
		{
			iIntegralSymmetry.col(i) = mVectorsVerticesBasisInv * rqiVectorsC[iVertices[iPermutation[iLinearlyIndependantVertices[i]]]];
		}
		
		Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > mCanBasis( mVectorsVerticesBasis * iIntegralSymmetry * mVectorsVerticesBasisInv );
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

bool QuadraticInteger_InfiniteNSymetries::FindIntegralSymmetryFromSubgraph( const vector< unsigned int >& iVertices )
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

bool QuadraticInteger_InfiniteNSymetries::bDottedSameWeight( const unsigned int& v1, const unsigned int& v2, const unsigned int& w1, const unsigned int& w2 ) const
{
	return ( qiDottedWeights[v1][v2] == qiDottedWeights[w1][w2] );
}

void QuadraticInteger_InfiniteNSymetries::WorkWithIsomorphisms( const vector< unsigned int >& iVertices, const vector< GraphInvolution >& iIsomorphisms )
{
	unsigned int i, j, k, iVerticesCount( iVertices.size() );
	
	// -------------------------------------------
	// Vectors corresponding to the vertices
	Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > m;
	m.resize( iVectorSize, iVerticesCount );
	
	for( unsigned int i(0); i < iVerticesCount; i++ )
		m.col(i) = rqiVectorsC[i];
	
	// --------------------------------------------
	// For each isomorphism
	for( auto iIso : iIsomorphisms )
	{
		// ------------------------------------------------
		// The column of m span the space of fixed points
		Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > m;
		m.resize( iVectorSize, iIso.iOrbitCount );
		
		j = 0;
		for( i = 0; i < iVerticesCount && j < iIso.iOrbitCount; i++ )
		{
			if( iIso.iPermutation[i] == i )
				m.col(j++) = rqiVectorsC[iVertices[i]];
			else if( iIso.iPermutation[i] < i )
				m.col(j++) = rqiVectorsC[iVertices[i]] + rqiVectorsC[ iVertices[iIso.iPermutation[i]]];
		}
		
		// We get a basis of the space of fixed points
		FullPivLU<Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic >> lu(m);
		Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > basisFixedPoints( lu.image(m) );
		
		unsigned int iFixedPointDimTemp( basisFixedPoints.cols() );
		if( iFixedPointDimTemp == iVectorSize ) // Soooo much fixed points
			continue;
		
		#pragma omp critical
		{
			if( !rqiBasisFixedPoints.cols() ) // These are the first fixed points
			{
				rqiBasisFixedPoints = basisFixedPoints;
				iFixedPointsDimension = iFixedPointDimTemp;
				usefulInvolutions.push_back( iIso );
			}
			else // We already have fixed points
			{
				// The goal here is to take the intersection of the previsous fixed points space and the new one
				Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > tempVectors( rqiBasisFixedPoints );
				
				tempVectors.conservativeResize( iVectorSize, tempVectors.cols() + iFixedPointDimTemp );
				for( j = 0; j < iFixedPointDimTemp; j++ )
					tempVectors.col(j + iFixedPointsDimension) = basisFixedPoints.col(j);
				
				FullPivLU<Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > > luTemp(tempVectors);
				Matrix< Rational<QuadraticIntegerBig>, Dynamic, Dynamic > ker( luTemp.kernel() );
				
				if( luTemp.dimensionOfKernel() == 0 )
				{
					bFinished = true;
					iFixedPointsDimension = 0;
					rqiBasisFixedPoints.resize( iVectorSize, 0 );
					usefulInvolutions.push_back( iIso );
				}
				
				if( ker.cols() < iFixedPointsDimension && !bFinished ) // If we earned something
				{
					usefulInvolutions.push_back( iIso );
					unsigned int iIntersectionDimension( ker.cols() );
					rqiBasisFixedPoints.resize( iVectorSize, iIntersectionDimension );
									
					for( j = 0; j < iIntersectionDimension; j++ )
					{
						Matrix< Rational<QuadraticIntegerBig>, Dynamic, 1 > v( ker(0,j)*tempVectors.col(0) );
						
						for( k = 1; k < iFixedPointsDimension; k++ )
							v += ker(k,j)*tempVectors.col(k);
						
						rqiBasisFixedPoints.col(j) = v;
					}
					
					iFixedPointsDimension = ker.cols();
					
					if( iFixedPointsDimension == 1 )
					{	
						if( rqiVectorNorm( rqiBasisFixedPoints.col(0) ) >= 0 )
						{
							bFinished = true;
						}
					}
					else if( iFixedPointsDimension == 2 )
					{
						// If the basis has two orthogonal vector of positive norm
						if( rqiVectorNorm( rqiBasisFixedPoints.col(0) ) >= 0 && rqiVectorNorm( rqiBasisFixedPoints.col(1) ) >= 0 && rqiVectorsProduct( rqiBasisFixedPoints.col(0), rqiBasisFixedPoints.col(1) ) == 0 )
							bFinished = true;
					}
				}
			}
		}
		
		if( bFinished )
			return;
	}
}

void QuadraticInteger_InfiniteNSymetries::print_basisFixedPoints( const string& strSpacer ) const
{
	for( unsigned int i(0); i < iFixedPointsDimension; i++ )
		cout << strSpacer << Matrix<Rational<QuadraticIntegerBig>, 1, Eigen::Dynamic>( rqiBasisFixedPoints.col(i).transpose() ) << endl;
}
