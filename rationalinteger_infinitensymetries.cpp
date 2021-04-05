#include "rationalinteger_infinitensymetries.h"

RationalInteger_InfiniteNSymetries::RationalInteger_InfiniteNSymetries(
    AlVin *alvin)
    : InfiniteNSymetries(alvin) {
  unsigned int i, j;

  auto vectors(alvin->get_vectors());

  riVectorsProducts = vector<vector<RationalInteger>>(
      vectorsCount, vector<RationalInteger>(vectorsCount, 0));

  // -------------------------------------------------
  // Quadratic form
  for (auto ai : qf)
    rriQF.push_back(
        Rational<RationalInteger>(*dynamic_cast<RationalInteger *>(ai)));

  // -------------------------------------------------
  // Vectors
  for (i = 0; i < vectorsCount; i++) {
    auto v(vector<Rational<RationalInteger>>(0));

    Matrix<Rational<RationalInteger>, Dynamic, 1> v2;

    for (j = 0; j <= dimension; j++)
      v.push_back(Rational<RationalInteger>(
          *dynamic_cast<RationalInteger *>(vectors[i][j])));

    v2 = Matrix<Rational<RationalInteger>, 1, Dynamic>::Map(&v[0], 1,
                                                            dimension + 1);
    rriVectorsC.push_back(v2);
  }

  // -------------------------------------------------
  // Products
  for (i = 0; i < vectorsCount; i++) {
    AlgebraicInteger *norm2(alvin->bilinearProduct(vectors[i], vectors[i]));
    riVectorsProducts[i][i] = *dynamic_cast<RationalInteger *>(norm2);
    delete norm2;

    for (j = i + 1; j < vectorsCount; j++) {
      AlgebraicInteger *product(alvin->bilinearProduct(vectors[i], vectors[j]));
      riVectorsProducts[i][j] = riVectorsProducts[j][i] =
          *dynamic_cast<RationalInteger *>(product);
      delete product;
    }
  }

  // -------------------------------------------------
  // Weights of dotted lines
  riDottedWeights = vector<vector<Rational<RationalInteger>>>(
      vectorsCount, vector<Rational<RationalInteger>>(vectorsCount, 0));
  for (i = 0; i < vectorsCount; i++) {
    for (j = i + 1; j < vectorsCount; j++) {
      if (graphMatrix[i][j] == 2) // dotted
      {
        auto num(alvin->bilinearProduct(vectors[i], vectors[j]));
        auto norm1(alvin->bilinearProduct(vectors[i], vectors[i]));
        auto norm2(alvin->bilinearProduct(vectors[j], vectors[j]));

        num->multiplyBy(num);
        norm1->multiplyBy(norm2);

        riDottedWeights[i][j] = riDottedWeights[j][i] =
            Rational<RationalInteger>(*dynamic_cast<RationalInteger *>(num),
                                      *dynamic_cast<RationalInteger *>(norm1));

        delete num;
        delete norm1;
        delete norm2;
      }
    }
  }
}

RationalInteger_InfiniteNSymetries::~RationalInteger_InfiniteNSymetries() {}

Rational<RationalInteger> RationalInteger_InfiniteNSymetries::rriVectorNorm(
    const Matrix<Rational<RationalInteger>, -1, 1> &rriV) const {
  if (rriV.rows() != iVectorSize)
    throw(string(
        "RationalInteger_InfiniteNSymetries::rriVectorNorm: Invalid size"));

  Rational<RationalInteger> rriRes(-rriQF[0] * rriV(0, 0) * rriV(0, 0));

  for (unsigned int j(1); j < iVectorSize; j++)
    rriRes += rriQF[j] * rriV(j, 0) * rriV(j, 0);

  return rriRes;
}

Rational<RationalInteger> RationalInteger_InfiniteNSymetries::rriVectorsProduct(
    const Matrix<Rational<RationalInteger>, Dynamic, 1> &rriV1,
    const Matrix<Rational<RationalInteger>, Dynamic, 1> &rriV2) const {
  Rational<RationalInteger> rriRes(-rriQF[0] * rriV1(0, 0) * rriV2(0, 0));

  for (unsigned int j(1); j < iVectorSize; j++)
    rriRes += rriQF[j] * rriV1(j, 0) * rriV2(j, 0);

  return rriRes;
}

vector<GraphInvolution>
RationalInteger_InfiniteNSymetries::FindIsomorphismsInSubgraph(
    const vector<unsigned int> &iVertices) {
  vector<GraphInvolution> iIsomorphisms;

  unsigned int iVerticesCount(iVertices.size()), i, j, k;

  // -------------------------------------------------------------
  // Do the corresponding vectors span the whole space?
  Matrix<Rational<RationalInteger>, Dynamic, Dynamic> mVectorsVertices_temp;
  mVectorsVertices_temp.resize(iVectorSize, iVerticesCount);

  for (i = 0; i < iVerticesCount; i++)
    mVectorsVertices_temp.col(i) = rriVectorsC[iVertices[i]];

  FullPivLU<Matrix<Rational<RationalInteger>, Dynamic, Dynamic>>
      lu_mVectorsVertices_temp(mVectorsVertices_temp);
  if (lu_mVectorsVertices_temp.rank() < iVectorSize)
    return iIsomorphisms;

  // -------------------------------------------------------------
  // We extract a basis of the space, consisting of a subset of vectors
  // designated by iVertices
  auto iPerm(lu_mVectorsVertices_temp.permutationQ().indices());

  vector<unsigned int> iLinearlyIndependantVertices(
      iVectorSize,
      0); ///< Strictly speaking, the indices i (in iVertices!) such that the
          ///< set vectors[ iVertices[i] ] is linearly indpendant

  std::sort(iPerm.derived().data(), iPerm.derived().data() + iVectorSize);
  Matrix<Rational<RationalInteger>, Dynamic, Dynamic> mVectorsVerticesBasis;
  mVectorsVerticesBasis.resize(iVectorSize, iVectorSize);
  for (i = 0; i < iVectorSize; i++) {
    iLinearlyIndependantVertices[i] = iPerm[i];
    mVectorsVerticesBasis.col(i) = mVectorsVertices_temp.col(iPerm[i]);
  }

  // -------------------------------------------------------------
  // Constructing the graph
  igraph_t graph;
  igraph_matrix_t adjMat;
  igraph_matrix_init(&adjMat, iVerticesCount, iVerticesCount);

  // vertices coloring
  igraph_vector_int_t verticesColors;
  igraph_vector_int_init(&verticesColors, iVerticesCount);

  for (i = 0; i < iVerticesCount; i++) {
    VECTOR(verticesColors)
    [i] = riVectorsProducts[iVertices[i]][iVertices[i]].val;

    for (j = i + 1; j < iVerticesCount; j++)
      MATRIX(adjMat, i, j) = graphMatrix[iVertices[i]][iVertices[j]];
  }

  if (igraph_adjacency(&graph, &adjMat, IGRAPH_ADJ_UNDIRECTED) !=
      IGRAPH_SUCCESS)
    throw(string(
        "InfiniteNSymetries::FindInvolutionInSubgraph: Error creating graph"));

  // -------------------------------------------------------------
  // Finding isomorphisms
  igraph_vector_ptr_t maps;
  igraph_vector_ptr_init(&maps, 1);
  igraph_get_isomorphisms_vf2(&graph, &graph, &verticesColors, &verticesColors,
                              0, 0, &maps, 0, 0, 0);
  long int iIsomorphismsCount(igraph_vector_ptr_size(&maps));

  for (int is(1); is < iIsomorphismsCount; is++) // Don't consider the identity
  {
    bool bIsInvolution(true), bIsIsomorphism(true);
    vector<unsigned int> iPermutation(iVerticesCount, -1);
    unsigned int iFixedPoints(0);

    igraph_vector_t *v((igraph_vector_t *)igraph_vector_ptr_e(&maps, is));

    // ---------------------------------------------------
    // Checking that it is an involution and getting the permutation
    for (k = 0; k < iVerticesCount; k++) {
      iPermutation[k] = (unsigned int)VECTOR(*v)[k];

      if (iPermutation[k] == k)
        iFixedPoints++;

      if ((unsigned int)VECTOR(*v)[(unsigned int)VECTOR(*v)[k]] != k) {
        bIsInvolution = false;
        break;
      }
    }

    igraph_vector_destroy(v);
    free(v);

    if (!bIsInvolution)
      continue;

    // ---------------------------------------------------
    // Checking that it is an isomorphism
    for (i = 0; i < iVerticesCount && bIsIsomorphism; i++) {
      for (j = i + 1; j < iVerticesCount; j++) {
        // Weight not preserved or dotted and weight of the dotted not preserved
        if ((graphMatrix[iVertices[i]][iVertices[j]] !=
             graphMatrix[iVertices[iPermutation[i]]]
                        [iVertices[iPermutation[j]]]) ||
            (graphMatrix[iVertices[i]][iVertices[j]] == 2 &&
             !bDottedSameWeight(iVertices[i], iVertices[j],
                                iVertices[iPermutation[i]],
                                iVertices[iPermutation[j]]))) {
          bIsIsomorphism = false;
          break;
        }

        // If does not preserve some angle
        if (riVectorsProducts[iVertices[i]][iVertices[j]] !=
            riVectorsProducts[iVertices[iPermutation[i]]]
                             [iVertices[iPermutation[j]]]) {
          bIsIsomorphism = false;
          break;
        }
      }
    }

    if (!bIsIsomorphism) // Does not preserve the dotted
      continue;

    // ---------------------------------------------------
    // Here, we still have to check it preserves the lattice generated by the
    // canonical basis
    Matrix<Rational<RationalInteger>, Dynamic, Dynamic>
        iIntegralSymmetry; ///< Permutation matrix wrt to a basis contained in
                           ///< the vectors corresponding to iVertices
    Matrix<Rational<RationalInteger>, Dynamic, Dynamic>
        mVectorsVerticesBasisInv(mVectorsVerticesBasis.inverse());

    iIntegralSymmetry.resize(iVectorSize, iVectorSize);
    for (i = 0; i < iVectorSize; i++)
      iIntegralSymmetry.col(i) =
          mVectorsVerticesBasisInv *
          rriVectorsC[iVertices[iPermutation[iLinearlyIndependantVertices[i]]]];

    Matrix<Rational<RationalInteger>, Dynamic, Dynamic> mCanBasis(
        mVectorsVerticesBasis * iIntegralSymmetry * mVectorsVerticesBasisInv);
    bool bPreservesLattice(true);

    for (i = 0; i < iVectorSize && bPreservesLattice; i++) {
      for (j = 0; j < iVectorSize; j++) {
        if (!mCanBasis(i, j).get_hasDenominatorOne()) {
          bPreservesLattice = false;
          break;
        }
      }
    }

    if (bPreservesLattice)
      iIsomorphisms.push_back({iVertices, iPermutation,
                               iLinearlyIndependantVertices,
                               (iVerticesCount + iFixedPoints) / 2});
  }

  // ----------------------------------------
  // Freeing memory
  igraph_destroy(&graph);
  igraph_matrix_destroy(&adjMat);

  igraph_vector_destroy(
      (igraph_vector_t *)igraph_vector_ptr_e(&maps, 0));  // identity
  free((igraph_vector_t *)igraph_vector_ptr_e(&maps, 0)); // identity

  igraph_vector_ptr_destroy(&maps);
  igraph_vector_int_destroy(&verticesColors);

  return iIsomorphisms;
}

bool RationalInteger_InfiniteNSymetries::FindIntegralSymmetryFromSubgraph(
    const vector<unsigned int> &iVertices) {
  unsigned int i, j;
  unsigned int iVerticesCount(iVertices.size());

  // -------------------------------------------------------------
  // Does the corresponding polyhedron have at least one vertex?
  vector<vector<unsigned int>> iCox(vector<vector<unsigned int>>(
      iVerticesCount, vector<unsigned int>(iVerticesCount, 2)));
  for (i = 0; i < iVerticesCount; i++) {
    for (j = i + 1; j < iVerticesCount; j++)
      iCox[i][j] = iCox[j][i] = coxeterMatrix[iVertices[i]][iVertices[j]];
  }

  CoxIter ci(iCox, dimension);
  ci.exploreGraph();
  ci.computeGraphsProducts();
  ci.computeEulerCharacteristicFVector();
  vector<unsigned int> iFV(ci.get_fVector());
  if (!iFV[0])
    return isFinished;

  // -------------------------------------------------------------
  // Finding the isomorphisms
  vector<GraphInvolution> iIsomorphisms(FindIsomorphismsInSubgraph(iVertices));
  if (!iIsomorphisms.size())
    return isFinished;

  // -------------------------------------------------------------
  // Going through
  WorkWithIsomorphisms(iVertices, iIsomorphisms);

  return isFinished;
}

bool RationalInteger_InfiniteNSymetries::bDottedSameWeight(
    const unsigned int &v1, const unsigned int &v2, const unsigned int &w1,
    const unsigned int &w2) const {
  return (riDottedWeights[v1][v2] == riDottedWeights[w1][w2]);
}

void RationalInteger_InfiniteNSymetries::WorkWithIsomorphisms(
    const vector<unsigned int> &iVertices,
    const vector<GraphInvolution> &iIsomorphisms) {
  unsigned int i, j, k, iVerticesCount(iVertices.size());

  // -------------------------------------------
  // Vectors corresponding to the vertices
  Matrix<Rational<RationalInteger>, Dynamic, Dynamic> m;
  m.resize(iVectorSize, iVerticesCount);

  for (unsigned int i(0); i < iVerticesCount; i++)
    m.col(i) = rriVectorsC[i];

  // --------------------------------------------
  // For each isomorphism
  for (auto iIso : iIsomorphisms) {
    // ------------------------------------------------
    // The column of m span the space of fixed points
    Matrix<Rational<RationalInteger>, Dynamic, Dynamic> m;
    m.resize(iVectorSize, iIso.iOrbitCount);

    j = 0;
    for (i = 0; i < iVerticesCount && j < iIso.iOrbitCount; i++) {
      if (iIso.iPermutation[i] == i)
        m.col(j++) = rriVectorsC[iVertices[i]];
      else if (iIso.iPermutation[i] < i)
        m.col(j++) = rriVectorsC[iVertices[i]] +
                     rriVectorsC[iVertices[iIso.iPermutation[i]]];
    }

    // We get a basis of the space of fixed points
    FullPivLU<Matrix<Rational<RationalInteger>, Dynamic, Dynamic>> lu(m);
    Matrix<Rational<RationalInteger>, Dynamic, Dynamic> basisFixedPoints(
        lu.image(m));

    unsigned int iFixedPointDimTemp(basisFixedPoints.cols());
    if (iFixedPointDimTemp == iVectorSize) // Soooo much fixed points
      continue;

#pragma omp critical
    {
      if (!rriBasisFixedPoints.cols()) // These are the first fixed points
      {
        rriBasisFixedPoints = basisFixedPoints;
        fixedPointsDimension = iFixedPointDimTemp;
        usefulInvolutions.push_back(iIso);
      } else // We already have fixed points
      {
        // The goal here is to take the intersection of the previsous fixed
        // points space and the new one
        Matrix<Rational<RationalInteger>, Dynamic, Dynamic> tempVectors(
            rriBasisFixedPoints);

        tempVectors.conservativeResize(iVectorSize,
                                       tempVectors.cols() + iFixedPointDimTemp);
        for (j = 0; j < iFixedPointDimTemp; j++)
          tempVectors.col(j + fixedPointsDimension) = basisFixedPoints.col(j);

        FullPivLU<Matrix<Rational<RationalInteger>, Dynamic, Dynamic>> luTemp(
            tempVectors);
        Matrix<Rational<RationalInteger>, Dynamic, Dynamic> ker(
            luTemp.kernel());

        if (luTemp.dimensionOfKernel() == 0) {
          isFinished = true;
          fixedPointsDimension = 0;
          rriBasisFixedPoints.resize(iVectorSize, 0);
          usefulInvolutions.push_back(iIso);
        }

        if (ker.cols() < fixedPointsDimension &&
            !isFinished) // If we earned something
        {
          usefulInvolutions.push_back(iIso);
          unsigned int iIntersectionDimension(ker.cols());
          rriBasisFixedPoints.resize(iVectorSize, iIntersectionDimension);

          for (j = 0; j < iIntersectionDimension; j++) {
            Matrix<Rational<RationalInteger>, Dynamic, 1> v(ker(0, j) *
                                                            tempVectors.col(0));

            for (k = 1; k < fixedPointsDimension; k++)
              v += ker(k, j) * tempVectors.col(k);

            rriBasisFixedPoints.col(j) = v;
          }

          fixedPointsDimension = ker.cols();

          if (fixedPointsDimension == 1) {
            if (rriVectorNorm(rriBasisFixedPoints.col(0)) >= 0)
              isFinished = true;
          } else if (fixedPointsDimension == 2) {
            // If the basis has two orthogonal vector of positive norm
            if (rriVectorNorm(rriBasisFixedPoints.col(0)) >= 0 &&
                rriVectorNorm(rriBasisFixedPoints.col(1)) >= 0 &&
                rriVectorsProduct(rriBasisFixedPoints.col(0),
                                  rriBasisFixedPoints.col(1)) == 0)
              isFinished = true;
          }
        }
      }
    }

    if (isFinished)
      return;
  }
}

void RationalInteger_InfiniteNSymetries::print_basisFixedPoints(
    const string &strSpacer) const {
  if (fixedPointsDimension == iVectorSize) {
    cout << strSpacer << "The fixed point space is full dimensional" << endl;
    return;
  }

  for (unsigned int i(0); i < fixedPointsDimension; i++)
    cout << strSpacer
         << Matrix<Rational<RationalInteger>, 1, Eigen::Dynamic>(
                rriBasisFixedPoints.col(i).transpose())
         << endl;
}
