#include "rationalinteger_alvinfractions.h"

RationalInteger_VFs::RationalInteger_VFs(
    vector<AlgebraicInteger *> possibleNorms2)
    : AlVinFractions(possibleNorms2),
      possibleNorms2Max(dynamic_cast<RationalInteger *>(
                            possibleNorms2[possibleNorms2.size() - 1])
                            ->get_iValue()) {}

RationalInteger_VFs::~RationalInteger_VFs() {}

void RationalInteger_VFs::computeNextAlVinFractions() {
  unsigned int iMin(lastMaximum);
  unsigned int iNumerator;

  lastMaximum = (iMin + batchSize) * (iMin + batchSize);

  // We generate fractions f such that: iMin < f <= iLastMaximum

  for (const auto &norm2 : possibleNorms2) {
    unsigned int iE(dynamic_cast<RationalInteger *>(norm2)->get_iValue());

    // ceil( isqrt(x) ) = isqrt(x - 1) + 1
    unsigned int xMin(iMin ? integerSqrt(iE * iMin - 1) + 1 : 1);
    unsigned int xMax(integerSqrt(iE * lastMaximum));

    if (xMin * xMin == iE * iMin) // because we want iMin < f and not iMin <= f
      xMin++;

    for (; xMin <= xMax; xMin++) {
      iNumerator = xMin * xMin * possibleNorms2Max / iE;

      AlVinFraction *vf =
          new AlVinFraction(new RationalInteger(xMin), new RationalInteger(iE),
                            new RationalInteger(iNumerator));

      alvinFractions.insert(lower_bound(alvinFractions.begin(),
                                        alvinFractions.end(), vf,
                                        isLessThanPtrAlVinFraction),
                            vf);
    }
  }

  alvinfractions_it = alvinFractions.begin();
}
