#include "infinitensymetries.h"

InfiniteNSymetries::InfiniteNSymetries(AlVin *alvin)
    : alvin(alvin), dimension(alvin->get_dimension()),
      iVectorSize(alvin->get_dimension() + 1),
      coxeterMatrix(alvin->get_coxeterMatrix()), isFinished(false),
      qf(alvin->get_qf()), fixedPointsDimension(alvin->get_dimension() + 1) {
  vectorsCount = coxeterMatrix.size();
  graphMatrix = vector<vector<unsigned int>>(
      vectorsCount, vector<unsigned int>(vectorsCount, 0));

  for (unsigned int i(0); i < vectorsCount; i++) {
    for (unsigned int j(i + 1); j < vectorsCount; j++)
      graphMatrix[i][j] = graphMatrix[j][i] =
          coxeterMatrix[i][j] >= 3
              ? coxeterMatrix[i][j]
              : (coxeterMatrix[i][j] == 2 ? 0 : coxeterMatrix[i][j] + 1);
  }
}

bool InfiniteNSymetries::Run(const unsigned int &iNRMin,
                             const unsigned int &iNRMax) {
  string bitmask;
  unsigned int i, j, N(vectorsCount);

  isFinished = false;

#pragma omp parallel for schedule(dynamic, 1)                                  \
    shared(N, iNRMax, iNRMin) private(bitmask, i, j) collapse(2)
  for (unsigned int K = iNRMin; K <= iNRMax; K++) {
    for (unsigned int iFirstComponent = 0; iFirstComponent <= N;
         ++iFirstComponent) {
      if (iFirstComponent > N - K)
        continue;

      vector<unsigned int> iSubset(K, 0);
      bitmask = std::string(iFirstComponent, 0); // iFirstComponent-1 0's
      bitmask.resize(iFirstComponent + K, 1);    // K 1's
      bitmask.resize(N, 0);                      // N-K trailing 0's

      do {
        if (isFinished)
          break;

        j = 0;
        for (i = 0; i < N; ++i) {
          if (bitmask[i])
            iSubset[j++] = i;
        }

        try {
          if (!isFinished && FindIntegralSymmetryFromSubgraph(iSubset)) {
            isFinished = true;
            break;
          }
        } catch (string strE) {
          cout << "-----------------------------\nEXCEPTION:\n" << strE << endl;
          cout << "Set: " << implode(",", iSubset) << endl;
          exit(0);
        }

      } while (std::prev_permutation(bitmask.begin(), bitmask.end()) &&
               bitmask[iFirstComponent]);
    }
  }

  return isFinished;
}

vector<GraphInvolution> InfiniteNSymetries::get_usefulInvolutions() const {
  return usefulInvolutions;
}

unsigned int InfiniteNSymetries::get_iFixedPointsDimension() const {
  return fixedPointsDimension;
}

InfiniteNSymetries::~InfiniteNSymetries() {}
