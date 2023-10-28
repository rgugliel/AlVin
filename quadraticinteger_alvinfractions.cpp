#include "quadraticinteger_alvinfractions.h"

QuadraticInteger_VFs::QuadraticInteger_VFs(
    vector<AlgebraicInteger *> aiPossibleNorms2,
    const QuadraticInteger &qiAlpha0)
    : qiAlpha0(qiAlpha0), AlVinFractions(aiPossibleNorms2),
      qiPossibleNorms2_max(dynamic_cast<QuadraticInteger *>(
          aiPossibleNorms2[aiPossibleNorms2.size() - 1])) {
  iBatchSize = 1;
}

QuadraticInteger_VFs::~QuadraticInteger_VFs() {}

void QuadraticInteger_VFs::computeNextAlVinFractions() {
  unsigned int iMin(iLastMaximum),
      iSizeTemp((iMin + iBatchSize) * (iMin + iBatchSize) - iLastMaximum);

  iLastMaximum += min(iSizeTemp, (unsigned int)1000000);

  // We generate fractions f such that: iMin < f <= iLastMaximum
  for (auto aiE : aiPossibleNorms2) {
    QuadraticInteger *qiE(dynamic_cast<QuadraticInteger *>(aiE));
    long int y, iYMax, iYMin, x, iXMin, iXMax, iSqrtEp, iTemp, iSqrtY2d;

    QuadraticInteger qiMinEpsilon(*qiE), qiMaxEpsilon(*qiE),
        qiConjEpsilon(*qiE);
    qiMinEpsilon.multiplyBy(iMin);
    qiMaxEpsilon.multiplyBy(iLastMaximum);
    qiConjEpsilon.conjugate();

    if (QuadraticInteger::isOneMod4) {
      // iSqrtEp = isqrt( \bar epsilon / ( -d * \bar \alpha_0 ) )
      QuadraticInteger qiTemp(*qiE);
      qiTemp.conjugate();
      QuadraticInteger qiDen(qiAlpha0);
      qiDen.conjugate();
      qiDen.multiplyBy(-QuadraticInteger::d);
      iSqrtEp = QuadraticInteger::iSQRT_quotient(qiTemp, qiDen);

      // isqrt( iLastMaximum * epsilon / d )
      iYMax = iSqrtEp +
              QuadraticInteger::iSQRT_quotient(
                  qiMaxEpsilon, QuadraticInteger(QuadraticInteger::d)) +
              1;

      // isqrt( iMin * epsilon / d )
      iYMin = QuadraticInteger::iSQRT_quotient(
                  qiMinEpsilon, QuadraticInteger(QuadraticInteger::d)) -
              iSqrtEp;

      qiTemp.set(&qiMinEpsilon);
      qiTemp.multiplyBy(4);

      long int iSQRTi4MinEpsilon(
          QuadraticInteger::iSQRT_quotient(qiTemp, QuadraticInteger(1)));

      long int iSQRTMaxEpsilon(
          QuadraticInteger::iSQRT_quotient(qiMaxEpsilon, QuadraticInteger(1)));

      // TODO: est ce que y < ymax et x < xmax --> pas de vérifications à faire?

      for (y = iYMin; y <= iYMax; y++) {
        iSqrtY2d = integerSqrt((unsigned long int)(y * y * QuadraticInteger::d));
        ;
        iTemp = iSQRTi4MinEpsilon - y - iSqrtY2d - 1;

        iXMin = (iTemp % 2) ? iTemp / 2 + 1 : iTemp / 2;
        iXMax = iSQRTMaxEpsilon + 1 - (y + iSqrtY2d) / 2;

        for (x = iXMin; x <= iXMax; x++) {
          QuadraticInteger qi(QuadraticInteger(x, y)),
              qiSquare(QuadraticInteger(x, y));

          if (qi.isLessThan(0) || (x == 0 && y == 0))
            continue;

          qiSquare.multiplyBy(&qi);

          if (!qiMinEpsilon.isLessThan(qiSquare))
            continue;

          if (!qiSquare.isLessOEThan(qiMaxEpsilon))
            continue;

          qiSquare.multiplyBy(-1);
          qiSquare.multiplyBy(&qiAlpha0);
          qiSquare.conjugate();

          if (!qiSquare.isLessOEThan(qiConjEpsilon))
            continue;

          QuadraticInteger *qiNum(new QuadraticInteger(qi));
          qiNum->multiplyBy(&qi);
          qiNum->multiplyBy(qiPossibleNorms2_max);
          qiNum->divideBy(qiE);

          AlVinFraction *vf = new AlVinFraction(
              new QuadraticInteger(qi), new QuadraticInteger(*qiE), qiNum);
          alvinfractions.insert(lower_bound(alvinfractions.begin(),
                                            alvinfractions.end(), vf,
                                            isLessThanPtrAlVinFraction),
                                vf);
        }
      }
    } else {
      QuadraticInteger qiMinConjAlpha0(qiAlpha0);
      qiMinConjAlpha0.conjugate();
      qiMinConjAlpha0.multiplyBy(-1);

      long int iSQRTsupConjEpsilonMinConjAlpha0(
          QuadraticInteger::iSQRTsup_quotient(qiConjEpsilon, qiMinConjAlpha0));

      iTemp =
          QuadraticInteger::iSQRT_quotient(qiMinEpsilon, QuadraticInteger(1)) -
          iSQRTsupConjEpsilonMinConjAlpha0;
      iXMin = (iTemp % 2) ? iTemp / 2 + 1 : iTemp / 2; // ceil
      iXMax = (QuadraticInteger::iSQRTsup_quotient(qiMaxEpsilon,
                                                   QuadraticInteger(1)) +
               iSQRTsupConjEpsilonMinConjAlpha0) /
              2;

      for (x = iXMin; x <= iXMax; x++) {
        iTemp = sqrtQuotient((unsigned long int)x * x,
                              (unsigned long int)QuadraticInteger::d);
        iYMin = QuadraticInteger::iSQRT_quotient(
                    qiMinEpsilon, QuadraticInteger(QuadraticInteger::d)) -
                iTemp;
        iYMax = QuadraticInteger::iSQRT_quotient(qiMaxEpsilon,
                                                 QuadraticInteger::d) -
                iTemp;

        for (y = iYMin; y <= iYMax; y++) {
          QuadraticInteger qi(QuadraticInteger(x, y)),
              qiSquare(QuadraticInteger(x, y));

          if (qi.isLessThan(0) || (x == 0 && y == 0))
            continue;

          qiSquare.multiplyBy(&qi);

          if (!qiMinEpsilon.isLessThan(qiSquare))
            continue;

          if (!qiSquare.isLessOEThan(qiMaxEpsilon))
            continue;

          qiSquare.multiplyBy(-1);
          qiSquare.multiplyBy(&qiAlpha0);
          qiSquare.conjugate();

          if (!qiSquare.isLessOEThan(qiConjEpsilon))
            continue;

          QuadraticInteger *qiNum(new QuadraticInteger(qi));
          qiNum->multiplyBy(&qi);
          qiNum->multiplyBy(qiPossibleNorms2_max);
          qiNum->divideBy(qiE);

          AlVinFraction *vf = new AlVinFraction(
              new QuadraticInteger(qi), new QuadraticInteger(*qiE), qiNum);
          alvinfractions.insert(lower_bound(alvinfractions.begin(),
                                            alvinfractions.end(), vf,
                                            isLessThanPtrAlVinFraction),
                                vf);
        }
      }
    }
  }

  alvinfractions_it = alvinfractions.begin();
}
