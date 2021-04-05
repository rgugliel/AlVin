#include "quadraticinteger_alvinfractions.h"

QuadraticInteger_VFs::QuadraticInteger_VFs(
    vector<AlgebraicInteger *> possibleNorms2, const QuadraticInteger &alpha0)
    : alpha0(alpha0), AlVinFractions(possibleNorms2),
      possibleNorms2Max(dynamic_cast<QuadraticInteger *>(
          possibleNorms2[possibleNorms2.size() - 1])) {
  batchSize = 1;
}

QuadraticInteger_VFs::~QuadraticInteger_VFs() {}

void QuadraticInteger_VFs::computeNextAlVinFractions() {
  unsigned int iMin(lastMaximum),
      iSizeTemp((iMin + batchSize) * (iMin + batchSize) - lastMaximum);

  lastMaximum += min(iSizeTemp, (unsigned int)1000000);

  // We generate fractions f such that: iMin < f <= iLastMaximum
  for (const auto &norm2 : possibleNorms2) {
    QuadraticInteger *qiE(dynamic_cast<QuadraticInteger *>(norm2));
    long int y, yMax, yMin, x, xMin, xMax, sqrtEp, temp, sqrtY2d;

    QuadraticInteger minEpsilon(*qiE), maxEpsilon(*qiE), conjEpsilon(*qiE);
    minEpsilon.multiplyBy(iMin);
    maxEpsilon.multiplyBy(lastMaximum);
    conjEpsilon.conjugate();

    if (QuadraticInteger::bIsOneMod4) {
      // sqrtEp = isqrt( \bar epsilon / ( -d * \bar \alpha_0 ) )
      QuadraticInteger qiTemp(*qiE);
      qiTemp.conjugate();
      QuadraticInteger qiDen(alpha0);
      qiDen.conjugate();
      qiDen.multiplyBy(-QuadraticInteger::d);
      sqrtEp = QuadraticInteger::iSQRT_quotient(qiTemp, qiDen);

      // isqrt( iLastMaximum * epsilon / d )
      yMax = sqrtEp +
             QuadraticInteger::iSQRT_quotient(
                 maxEpsilon, QuadraticInteger(QuadraticInteger::d)) +
             1;

      // isqrt( iMin * epsilon / d )
      yMin = QuadraticInteger::iSQRT_quotient(
                 minEpsilon, QuadraticInteger(QuadraticInteger::d)) -
             sqrtEp;

      qiTemp.set(&minEpsilon);
      qiTemp.multiplyBy(4);

      long int iSQRTi4MinEpsilon(
          QuadraticInteger::iSQRT_quotient(qiTemp, QuadraticInteger(1)));

      long int iSQRTMaxEpsilon(
          QuadraticInteger::iSQRT_quotient(maxEpsilon, QuadraticInteger(1)));

      // TODO: est ce que y < ymax et x < xmax --> pas de vérifications à faire?

      for (y = yMin; y <= yMax; y++) {
        sqrtY2d = integerSqrt((unsigned long int)(y * y * QuadraticInteger::d));
        temp = iSQRTi4MinEpsilon - y - sqrtY2d - 1;

        xMin = (temp % 2) ? temp / 2 + 1 : temp / 2;
        xMax = iSQRTMaxEpsilon + 1 - (y + sqrtY2d) / 2;

        for (x = xMin; x <= xMax; x++) {
          QuadraticInteger qi(QuadraticInteger(x, y)),
              qiSquare(QuadraticInteger(x, y));

          if (qi.isLessThan(0) || (x == 0 && y == 0))
            continue;

          qiSquare.multiplyBy(&qi);

          if (!minEpsilon.isLessThan(qiSquare))
            continue;

          if (!qiSquare.isLessOEThan(maxEpsilon))
            continue;

          qiSquare.multiplyBy(-1);
          qiSquare.multiplyBy(&alpha0);
          qiSquare.conjugate();

          if (!qiSquare.isLessOEThan(conjEpsilon))
            continue;

          QuadraticInteger *qiNum(new QuadraticInteger(qi));
          qiNum->multiplyBy(&qi);
          qiNum->multiplyBy(possibleNorms2Max);
          qiNum->divideBy(qiE);

          AlVinFraction *vf = new AlVinFraction(
              new QuadraticInteger(qi), new QuadraticInteger(*qiE), qiNum);
          alvinFractions.insert(lower_bound(alvinFractions.begin(),
                                            alvinFractions.end(), vf,
                                            isLessThanPtrAlVinFraction),
                                vf);
        }
      }
    } else {
      QuadraticInteger qiMinConjAlpha0(alpha0);
      qiMinConjAlpha0.conjugate();
      qiMinConjAlpha0.multiplyBy(-1);

      long int iSQRTsupConjEpsilonMinConjAlpha0(
          QuadraticInteger::iSQRTsup_quotient(conjEpsilon, qiMinConjAlpha0));

      temp = QuadraticInteger::iSQRT_quotient(minEpsilon, QuadraticInteger(1)) -
             iSQRTsupConjEpsilonMinConjAlpha0;
      xMin = (temp % 2) ? temp / 2 + 1 : temp / 2; // ceil
      xMax = (QuadraticInteger::iSQRTsup_quotient(maxEpsilon,
                                                  QuadraticInteger(1)) +
              iSQRTsupConjEpsilonMinConjAlpha0) /
             2;

      for (x = xMin; x <= xMax; x++) {
        temp = sqrtQuotient((unsigned long int)x * x,
                            (unsigned long int)QuadraticInteger::d);
        yMin = QuadraticInteger::iSQRT_quotient(
                   minEpsilon, QuadraticInteger(QuadraticInteger::d)) -
               temp;
        yMax =
            QuadraticInteger::iSQRT_quotient(maxEpsilon, QuadraticInteger::d) -
            temp;

        for (y = yMin; y <= yMax; y++) {
          QuadraticInteger qi(QuadraticInteger(x, y)),
              qiSquare(QuadraticInteger(x, y));

          if (qi.isLessThan(0) || (x == 0 && y == 0))
            continue;

          qiSquare.multiplyBy(&qi);

          if (!minEpsilon.isLessThan(qiSquare))
            continue;

          if (!qiSquare.isLessOEThan(maxEpsilon))
            continue;

          qiSquare.multiplyBy(-1);
          qiSquare.multiplyBy(&alpha0);
          qiSquare.conjugate();

          if (!qiSquare.isLessOEThan(conjEpsilon))
            continue;

          QuadraticInteger *qiNum(new QuadraticInteger(qi));
          qiNum->multiplyBy(&qi);
          qiNum->multiplyBy(possibleNorms2Max);
          qiNum->divideBy(qiE);

          AlVinFraction *vf = new AlVinFraction(
              new QuadraticInteger(qi), new QuadraticInteger(*qiE), qiNum);
          alvinFractions.insert(lower_bound(alvinFractions.begin(),
                                            alvinFractions.end(), vf,
                                            isLessThanPtrAlVinFraction),
                                vf);
        }
      }
    }
  }

  alvinfractions_it = alvinFractions.begin();
}
