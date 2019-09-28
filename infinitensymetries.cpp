#include "infinitensymetries.h"

InfiniteNSymetries::InfiniteNSymetries(AlVin *alvin)
    : alvin(alvin), iDimension(alvin->get_iDimension()),
      iVectorSize(alvin->get_iDimension() + 1),
      iCoxeterMatrix(alvin->get_iCoxeterMatrix()), bFinished(false),
      aiQF(alvin->get_aiQF()),
      iFixedPointsDimension(alvin->get_iDimension() + 1) {
  iVectorsCount = iCoxeterMatrix.size();
  iGraphMatrix = vector<vector<unsigned int>>(
      iVectorsCount, vector<unsigned int>(iVectorsCount, 0));

  for (unsigned int i(0); i < iVectorsCount; i++) {
    for (unsigned int j(i + 1); j < iVectorsCount; j++)
      iGraphMatrix[i][j] = iGraphMatrix[j][i] =
          iCoxeterMatrix[i][j] >= 3
              ? iCoxeterMatrix[i][j]
              : (iCoxeterMatrix[i][j] == 2 ? 0 : iCoxeterMatrix[i][j] + 1);
  }
}

bool InfiniteNSymetries::Run(const unsigned int &iNRMin,
                             const unsigned int &iNRMax) {
  string bitmask;
  unsigned int i, j, N(iVectorsCount);

  bFinished = false;

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
        if (bFinished)
          break;

        j = 0;
        for (i = 0; i < N; ++i) {
          if (bitmask[i])
            iSubset[j++] = i;
        }

        try {
          if (!bFinished && FindIntegralSymmetryFromSubgraph(iSubset)) {
            bFinished = true;
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

  return bFinished;
}

vector<GraphInvolution> InfiniteNSymetries::get_usefulInvolutions() const {
  return usefulInvolutions;
}

unsigned int InfiniteNSymetries::get_iFixedPointsDimension() const {
  return iFixedPointsDimension;
}

InfiniteNSymetries::~InfiniteNSymetries() {}
