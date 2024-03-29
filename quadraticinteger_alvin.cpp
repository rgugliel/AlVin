#include "quadraticinteger_alvin.h"

QuadraticInteger_AlVin::QuadraticInteger_AlVin(
    const vector<QuadraticInteger> &iQuadraticFormCoeffs,
    const string &strOuputMathematicalFormat, const bool &bWriteInfo,
    const bool &bDebug)
    : AlVin(strOuputMathematicalFormat, bWriteInfo, bDebug),
      d(QuadraticInteger::d), dISqrt(integerSqrt((unsigned int)QuadraticInteger::d)),
      isOneMod4(QuadraticInteger::isOneMod4) {
  bool bNegativeFound(false);

  if (QuadraticInteger::d == 0)
    throw(string("d is not set"));

  // -------------------------------------------------
  // Checking the coefficients
  for (auto i : iQuadraticFormCoeffs) {
    if (i.isLessThan(0)) {
      if (!bNegativeFound) {
        QuadraticInteger *qi(new QuadraticInteger(-i.a, -i.b));

        QuadraticInteger qiTemp(*qi);
        qiTemp.conjugate();

        if (qiTemp.isGreaterThan(0))
          throw(string("Quadratic form is not admissible"));

        aiQF.insert(aiQF.begin(), qi);

        bNegativeFound = true;
      } else
        throw(string("Quadratic form has incorrect signature"));
    } else if (i.isGreaterThan(0)) {
      QuadraticInteger *qi(new QuadraticInteger(i));

      QuadraticInteger qiTemp(*qi);
      qiTemp.conjugate();
      if (qiTemp.isLessThan(0))
        throw(string("Quadratic form is not admissible"));

      aiQF.insert(
          lower_bound(aiQF.begin() + (bNegativeFound ? 1 : 0), aiQF.end(), qi),
          qi);
    }
  }

  if (!bNegativeFound)
    throw(string("Quadratic form has incorrect signature"));

  // -------------------------------------------------
  // global work
  initializations();
  qiVectorCurrent = vector<QuadraticInteger>(iDimension + 1, 0);

  // -------------------------------------------------
  // local copy of the quadratic form (to avoid dynamic_cast)
  for (auto i : aiQF) {
    QuadraticInteger *qi(dynamic_cast<QuadraticInteger *>(i));
    qiQF.push_back(QuadraticInteger(*qi));
  }

  // -------------------------------------------------
  // Possible values for (e,e)
  findPossibleNorms2();

  // -------------------------------------------------
  // Displaying some information
  if (bWriteInfo)
    print_initialInformation();
}

void QuadraticInteger_AlVin::findPossibleNorms2() {
  vector<AlgebraicInteger *> aiPossibleNorms2(
      {new QuadraticInteger(1)});     ///< Possible values for (e,e)
  vector<QuadraticInteger> qiFactors; ///< Prime numbers and fundamental unit

  // ----------------------------------------------------------
  // Compute for ever coefficient a of the quadratic form the prime numbers
  // which divide a
  for (unsigned int i(0); i <= iDimension; i++) {
    vector<QuadraticInteger> qiTemp(qiQF[i].qiPrimeFactors());
    qiFactors.insert(qiFactors.end(), qiTemp.begin(), qiTemp.end());
  }

  // Removing duplicates
  sort(qiFactors.begin(), qiFactors.end());
  qiFactors = vector<QuadraticInteger>(
      qiFactors.begin(), unique(qiFactors.begin(), qiFactors.end()));

  // Prime factors of 2
  vector<QuadraticInteger> qiTemp(QuadraticInteger::qiFactorsRationalPrime(
      2, true)); // here we want the two factors of two (in case it splits)
  qiFactors.insert(qiFactors.end(), qiTemp.begin(), qiTemp.end());

  // Fundamental unit
  qiFactors.push_back(
      QuadraticInteger(QuadraticInteger::iFundamentalUnits[d][0],
                       QuadraticInteger::iFundamentalUnits[d][1]));

  unsigned int iFactorsCount(qiFactors.size());
  if (iFactorsCount > 8)
    throw(string("Number of values for (e,e) is too big"));

  // ----------------------------------------------------------
  // Powers of two
  vector<unsigned int> pow2({1});
  for (unsigned int i(1); i < iFactorsCount; i++)
    pow2[i] = pow2[i - 1] * 2;

  // ----------------------------------------------------------
  // Compute the products
  unsigned int iMax(pow(2, qiFactors.size()) - 1), j;

  for (unsigned int i = 1; i <= iMax; i++) {
    QuadraticInteger *qiProduct(new QuadraticInteger(1));

    for (j = 0; j < iFactorsCount; j++) {
      if (i & pow2[j])
        qiProduct->multiplyBy(&qiFactors[j]);
    }

    // is the product admissible (i.e. positive under galois embeddings)
    qiProduct->conjugate();
    if (qiProduct->isLessThan(0)) {
      qiProduct->conjugate();

      delete qiProduct;
      continue;
    }
    qiProduct->conjugate();

    // yes, we keep it
    aiPossibleNorms2.push_back(qiProduct);
  }

  // removing duplicates
  sort(aiPossibleNorms2.begin(), aiPossibleNorms2.end(),
       isLessThanPtrAlgebraicInteger);
  aiPossibleNorms2 = vector<AlgebraicInteger *>(
      aiPossibleNorms2.begin(),
      unique(aiPossibleNorms2.begin(), aiPossibleNorms2.end(),
             isEqualToPtrAlgebraicInteger));

  vf = new QuadraticInteger_VFs(aiPossibleNorms2, qiQF[0]);
}

bool QuadraticInteger_AlVin::PreRun() {
  if (d != QuadraticInteger::d)
    throw(string("QuadraticInteger_AlVin::The value of d has changed"));

  qiDQF.clear();
  qiDQFConj.clear();
  qi2QF.clear();
  for (unsigned int i(0); i <= iDimension; i++) {
    qiDQF.push_back(QuadraticInteger(qiQF[i]));
    qiDQF[i].multiplyBy(d);

    qiDQFConj.push_back(QuadraticInteger(qiDQF[i]));
    qiDQFConj[i].conjugate();

    qi2QF.push_back(qiQF[i]);
    qi2QF[i].multiplyBy(2);
  }

  return true;
}

void QuadraticInteger_AlVin::addVectorChild(
    const std::vector<AlgebraicInteger *> &aiVector) {
  vector<QuadraticInteger *> qiV;
  for (auto ai : aiVector)
    qiV.push_back(dynamic_cast<QuadraticInteger *>(ai));

  qiVectors.push_back(qiV);
}

int QuadraticInteger_AlVin::addVector_iFindWeight(
    AlgebraicInteger *aiNumerator, AlgebraicInteger *aiDenominator) {
  // -------------------------------------------
  // multiply both elements by the conjugate
  AlgebraicInteger *aiTemp(aiDenominator->copy());
  dynamic_cast<QuadraticInteger *>(aiTemp)->conjugate();

  aiDenominator->multiplyBy(aiTemp);
  aiNumerator->multiplyBy(aiTemp);

  delete aiTemp;

  // -------------------------------------------
  // simplification
  AlgebraicInteger *aiGCD(aiNumerator->copy());
  aiGCD->gcd(aiDenominator);

  aiNumerator->divideBy(aiGCD);
  aiDenominator->divideBy(aiGCD);

  if (aiNumerator->isLessThan(-1)) {
    aiNumerator->multiplyBy(-1);
    aiDenominator->multiplyBy(-1);
  }

  // -------------------------------------------
  // Proper tests
  if (aiNumerator->isEqualTo(1)) {
    if (QuadraticInteger::d == 2 &&
        aiDenominator->isEqualTo(QuadraticInteger(4, -2)))
      return 8;
    else if (QuadraticInteger::d == 5 &&
             aiDenominator->isEqualTo(QuadraticInteger(8, -4)))
      return 5;

    if (aiDenominator->isEqualTo(4))
      return 3;
    else if (aiDenominator->isEqualTo(2))
      return 4;
  }

  if (QuadraticInteger::d == 2) {
    if ((aiNumerator->isEqualTo(QuadraticInteger(1, 1)) &&
         aiDenominator->isEqualTo(QuadraticInteger(0, 2))))
      return 8;
  } else if (QuadraticInteger::d == 3) {
    if (aiNumerator->isEqualTo(QuadraticInteger(7, -4)) &&
        aiDenominator->isEqualTo(QuadraticInteger(28, -16)))
      return 3;
    else if (aiNumerator->isEqualTo(QuadraticInteger(2, -1)) &&
             aiDenominator->isEqualTo(QuadraticInteger(4, -2)))
      return 4;
    else if ((aiNumerator->isEqualTo(QuadraticInteger(2, -1)) &&
              aiDenominator->isEqualTo(QuadraticInteger(28, -16))) ||
             (aiNumerator->isEqualTo(1) &&
              aiDenominator->isEqualTo(QuadraticInteger(8, -4))))
      return 12;
  } else if (QuadraticInteger::d == 5) {
    if (aiNumerator->isEqualTo(QuadraticInteger(1, 1)) &&
        aiDenominator->isEqualTo(4))
      return 5;
    else if (aiNumerator->isEqualTo(QuadraticInteger(5, 5)) &&
             aiDenominator->isEqualTo(QuadraticInteger(8, 4)))
      return 10;
    else if (aiNumerator->isEqualTo(QuadraticInteger(2, 1)) &&
             aiDenominator->isEqualTo(4))
      return 10;
  } else if (QuadraticInteger::d == 7) {
    if (aiNumerator->isEqualTo(QuadraticInteger(127, -48)) &&
        aiDenominator->isEqualTo(QuadraticInteger(508, -192)))
      return 3;
    else if (aiNumerator->isEqualTo(QuadraticInteger(8, -3)) &&
             aiDenominator->isEqualTo(QuadraticInteger(16, -6)))
      return 4;
  }

  return -2;
}

void QuadraticInteger_AlVin::findVector(AlgebraicInteger *aiX0,
                                        AlgebraicInteger *aiNorm2) {
  QuadraticInteger *qi0(dynamic_cast<QuadraticInteger *>(aiX0));
  QuadraticInteger *qiNorm2(dynamic_cast<QuadraticInteger *>(aiNorm2));

  // -----------------------------------------------------
  // Preliminary work
  qiBilinearProducts = vector<QuadraticInteger>(iVectorsCount_second, 0);

  // initial value of iSumComp
  QuadraticInteger qiSumComp(*qiNorm2);
  QuadraticInteger qiTemp(*qi0);
  qiTemp.multiplyBy(qi0);
  qiTemp.multiplyBy(&qiQF[0]);
  qiSumComp.add(&qiTemp);

  // initial value of qiBilinearProducts
  for (unsigned int i(0); i < iVectorsCount_second; i++) {
    // = - qi[0] * k0 * qiVectors[i + iDimension][0]
    qiBilinearProducts[i].set(&qiQF[0]);
    qiBilinearProducts[i].multiplyBy(-1);
    qiBilinearProducts[i].multiplyBy(qi0);
    qiBilinearProducts[i].multiplyBy(qiVectors[i + iDimension][0]);
  }

  qiVectorCurrent[0].set(qi0);

  findVector(qi0, dynamic_cast<QuadraticInteger *>(aiNorm2), 1, qiSumComp,
             *qi0);
}

void QuadraticInteger_AlVin::findVector(QuadraticInteger *qi0,
                                        QuadraticInteger *qiNorm2,
                                        unsigned int iIndex,
                                        QuadraticInteger qiSumComp,
                                        QuadraticInteger qiGCDComponents) {
  if (qiSumComp.a == 0 && qiSumComp.b == 0) // We have a candidate
  {
    for (; iIndex <= iDimension; iIndex++)
      qiVectorCurrent[iIndex] = 0;

    if (qiGCDComponents.isInvertible())
      addCandidate();

    return;
  }

  /*
   * If d=1 (4), then iFirstMin=ymin, iFirstMax=ymax, iSecMin=xmin and
   * iSecMax=xmax. In particular, we first fix y and then x. If d=2,3 (4), then
   * iFirstMin=xmin, iFirstMax=xmax, iSecMin=ymin and iSecMax=ymax. In
   * particular, we first fix x and then y.
   */

  bool bAdmissible, bReturn(false);
  unsigned int j;
  long int iFirstMin, iFirstMax, iSec, iSecMin, iSecMax, iTemp, iTemp2;
  QuadraticInteger qi(0), qiSquare(0), qiTemp(0), qiPartialNorm(0),
      qiLastComponent(0), qiConjSumComp(qiSumComp);

  qiConjSumComp.conjugate();

  // -------------------------------------------------------------------
  // First bounds for the loop
  if (isOneMod4) {
    qiSumComp.conjugate();
    iFirstMax = QuadraticInteger::iSQRT_quotient(qiSumComp, qiDQFConj[iIndex]);

    iFirstMin = -iFirstMax - 1;

    qiSumComp.conjugate();
    iFirstMax += QuadraticInteger::iSQRT_quotient(qiSumComp, qiDQF[iIndex]) + 1;
  } else {
    qiSumComp.conjugate();
    iFirstMax = QuadraticInteger::iSQRT_quotient(qiSumComp, qi2QF[iIndex]);

    iFirstMin = -iFirstMax - 1;

    qiSumComp.conjugate();
    iFirstMax += QuadraticInteger::iSQRT_quotient(qiSumComp, qi2QF[iIndex]) + 1;
  }

  for (long int iFirst(iFirstMin); iFirst <= iFirstMax && !bReturn; iFirst++) {
    if (isOneMod4) {
      iTemp = iFirst + integerSqrt((unsigned long int)iFirst * iFirst * d) + 1;
      iSecMin = -iTemp / 2;

      iSecMax = QuadraticInteger::iSQRT_quotient(qiSumComp, qiQF[iIndex]) + 1 -
                ((iTemp % 2) ? iTemp / 2 + 1 : iTemp / 2);

      if (iComponentLessThan[iIndex]) // k_iIndex <= k_(iIndex-1) ?
      {
        iTemp2 = qiVectorCurrent[iComponentLessThan[iIndex]].b - iFirst;
        if (iTemp2 >= 0)
          iTemp =
              qiVectorCurrent[iComponentLessThan[iIndex]].a +
              (iTemp2 + sqrtSup<unsigned long int>(d * iTemp2 * iTemp2)) / 2;
        else {
          iTemp2 = abs(iTemp2);
          iTemp2 += integerSqrt<unsigned long int>(iTemp2 * iTemp2 * d);

          iTemp = qiVectorCurrent[iComponentLessThan[iIndex]].a -
                  ((iTemp2 % 2) ? iTemp2 / 2 + 1 : iTemp2 / 2);
        }

        iSecMax = min(iSecMax, iTemp);
      }
    } else {
      iTemp = sqrtQuotient((unsigned long int)iFirst * iFirst,
                            (unsigned long int)d);
      iSecMin = -iTemp - 1;
      iSecMax = QuadraticInteger::QuadraticInteger::iSQRT_quotient(
                    qiSumComp, qiDQF[iIndex]) -
                iTemp;

      if (iComponentLessThan[iIndex]) // k_iIndex <= k_(iIndex-1) ?
      {
        iTemp2 = qiVectorCurrent[iComponentLessThan[iIndex]].a - iFirst;

        if (iTemp2 >= 0)
          iTemp = sqrtQuotient<unsigned long int>(iTemp2 * iTemp2, d) +
                  qiVectorCurrent[iComponentLessThan[iIndex]].b;
        else
          iTemp = -sqrtSupQuotient<unsigned long int>(iTemp2 * iTemp2, d) +
                  qiVectorCurrent[iComponentLessThan[iIndex]].b;

        iSecMax = min(iSecMax, iTemp);
      }
    }

    for (iSec = iSecMin; iSec <= iSecMax; iSec++) {
      if (isOneMod4) {
        qi.a = iSec;
        qi.b = iFirst;
      } else {
        qi.a = iFirst;
        qi.b = iSec;
      }

      qiSquare.set(&qi);
      qiSquare.multiplyBy(&qi);

      // ------------------------------------------------------
      // crystallographic condition
      qiTemp.set(&qi2QF[iIndex]);
      qiTemp.multiplyBy(&qi);
      if (!qiTemp.isDivisibleBy(qiNorm2))
        continue;

      qiPartialNorm.set(&qiSquare);
      qiPartialNorm.multiplyBy(&qiQF[iIndex]);

      qiPartialNorm.conjugate();
      if (!qiPartialNorm.isLessOEThan(qiConjSumComp))
        continue;
      qiPartialNorm.conjugate();

      // ------------------------------------------------------
      // does iFirst + T * iSec satisfies the inequalities?
      if (!qiPartialNorm.isLessOEThan(qiSumComp))
        continue;

      if (!qi.isGreaterOEThan(0))
        continue;

      // TODO: incorporer la valeur MAX pour les composantes suite à
      // qiBilinearProducts
      if (iComponentLessThan[iIndex] &&
          !qi.isLessOEThan(qiVectorCurrent[iComponentLessThan[iIndex]]))
        continue;

      // ------------------------------------------------------
      // Is the vector admissible?
      bAdmissible = true;

      for (j = 0; j < iVectorsCount_second; j++) {
        // qiTemp = qiQF[iIndex] * ( qi - qiLastComponent ) * qiVectors[j +
        // iDimension][iIndex]
        qiTemp.set(&qi);
        qiTemp.substract(&qiLastComponent);
        qiTemp.multiplyBy(&qiQF[iIndex]);
        qiTemp.multiplyBy(qiVectors[j + iDimension][iIndex]);
        qiBilinearProducts[j].add(&qiTemp);

        // TODO: garder en mémoire la valeur max?

        if (qiBilinearProducts[j].isGreaterThan(
                0)) // qiVectorCurrent is not admissible
          bAdmissible = false;
      }

      qiLastComponent.set(&qi);

      if (!bAdmissible) {
        if (qi.a > 0 && qi.b > 0) // The products will grow
        {
          bReturn = true;
          break;
        }

        continue;
      }

      qiVectorCurrent[iIndex] = qi;

      if (iIndex == iDimension && qiSumComp.isEqualTo(qiPartialNorm)) {
        qiGCDComponents.gcd(&qi);
        if (qiGCDComponents.isInvertible())
          addCandidate();

        bReturn = true;
        break;
      }

      if (iIndex < iDimension) {
        QuadraticInteger qiSumSub(qiSumComp);
        qiSumSub.substract(&qiPartialNorm);

        qiTemp.set(&qiGCDComponents);
        qiTemp.gcd(&qi);

        findVector(qi0, qiNorm2, iIndex + 1, qiSumSub, qiTemp);
      }
    }
  }

  for (unsigned int i(0); i < iVectorsCount_second; i++) {
    qiTemp.set(&qiQF[iIndex]);
    qiTemp.multiplyBy(qiVectors[i + iDimension][iIndex]);
    qiTemp.multiplyBy(&qiLastComponent);
    qiBilinearProducts[i].substract(&qiTemp);
  }
}

void QuadraticInteger_AlVin::addCandidate() {
  vector<AlgebraicInteger *> aiV;
  for (auto i : qiVectorCurrent)
    aiV.push_back(new QuadraticInteger(i));

  aiVectors_candidates.push_back(aiV);
}

string QuadraticInteger_AlVin::get_strField() const {
  return string("Q[ sqrt(" + to_string(QuadraticInteger::d) + ") ]");
}
