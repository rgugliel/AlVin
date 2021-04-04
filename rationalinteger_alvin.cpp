#include "rationalinteger_alvin.h"

RationalInteger_AlVin::RationalInteger_AlVin(
    const vector<int> &iQuadraticFormCoeffs,
    const string &strOuputMathematicalFormat, const bool &bWriteInfo,
    const bool &bDebug)
    : AlVin(strOuputMathematicalFormat, bWriteInfo, bDebug) {
  bool bNegativeFound(false);

  // -------------------------------------------------
  // Checking the coefficients
  for (auto i : iQuadraticFormCoeffs) {
    if (i < 0 && !bNegativeFound) {
      RationalInteger *ri(new RationalInteger(-i));

      qf.insert(qf.begin(), ri);

      bNegativeFound = true;
    } else if (i < 0)
      throw(string("Quadratic form has incorrect signature"));
    else if (i > 0) {
      RationalInteger *ri(new RationalInteger(i));

      qf.insert(
          lower_bound(qf.begin() + (bNegativeFound ? 1 : 0), qf.end(), ri), ri);
    }
  }

  if (!bNegativeFound)
    throw(string("Quadratic form has incorrect signature"));

  // -------------------------------------------------
  // global work
  initializations();
  iVectorCurrent = vector<int>(dimension + 1, 0);

  // -------------------------------------------------
  // local copy of the quadratic form (computations performs faster without the
  // class-interface RationalInteger)
  for (auto i : qf) {
    RationalInteger *ri(dynamic_cast<RationalInteger *>(i));
    riQF.push_back(ri->get_iValue());
  }

  bLastTwoCoefficientsQFAre1Mod4 =
      riQF[dimension - 1] % 4 == 1 && riQF[dimension] % 4 == 1;

  // -------------------------------------------------
  // GCDs
  iGCDs = vector<unsigned int>(dimension + 1, 0);

  iGCDs[dimension - 1] = ugcd(riQF[dimension - 1], riQF[dimension]);
  if (iGCDs[dimension - 1] == 1)
    iGCDs[dimension - 1] = 0;

  for (unsigned int i(dimension - 2); i > 0; i--) {
    iGCDs[i] = ugcd(riQF[i], iGCDs[i + 1]);
    if (iGCDs[i] == 1) {
      iGCDs[i] = 0;
      break;
    }
  }

  // -------------------------------------------------
  // Possible values for (e,e)
  findPossibleNorms2();

  // -------------------------------------------------
  // Displaying some information
  if (bWriteInfo)
    print_initialInformation();
}

bool RationalInteger_AlVin::PreRun() { return true; }

void RationalInteger_AlVin::findPossibleNorms2() {
  vector<AlgebraicInteger *> aiPossibleNorms2(
      {new RationalInteger(1)}); ///< Possible values for (e,e)

  vector<unsigned int> iPrimeNumbers;

  // ----------------------------------------------------------
  // Compute for ever coefficient a of the quadratic form the prime numbers
  // which divide a
  for (unsigned int i(0); i <= dimension; i++) {
    unsigned int iDivisor(3), n(riQF[i]);
    if (n == 1)
      continue;

    if (!(n % 2)) {
      iPrimeNumbers.push_back(2);
      do {
        n /= 2;
      } while (!(n % 2));
    }

    while (n > 1) {
      if (!(n % iDivisor)) {
        iPrimeNumbers.push_back(iDivisor);
        do {
          n /= iDivisor;
        } while (!(n % iDivisor));
      }

      iDivisor += 2;
    }
  }

  sort(iPrimeNumbers.begin(), iPrimeNumbers.end());
  iPrimeNumbers =
      vector<unsigned int>(iPrimeNumbers.begin(),
                           unique(iPrimeNumbers.begin(), iPrimeNumbers.end()));

  iPrimeNumbers.push_back(2);

  unsigned int iPrimeNumbersCount(iPrimeNumbers.size());
  if (iPrimeNumbersCount > 8)
    throw(string("Number of values for (e,e) is too big"));

  // ----------------------------------------------------------
  // Powers of two
  vector<unsigned int> pow2({1});
  for (unsigned int i(1); i < iPrimeNumbersCount; i++)
    pow2[i] = pow2[i - 1] * 2;

  // ----------------------------------------------------------
  // Compute the products of prime numbers
  unsigned int iMax(pow(2, iPrimeNumbers.size()) - 1), j, iProduct;

  for (unsigned int i = 1; i <= iMax; i++) {
    iProduct = 1;
    for (j = 0; j < iPrimeNumbersCount; j++) {
      if (i & pow2[j])
        iProduct *= iPrimeNumbers[j];
    }

    aiPossibleNorms2.push_back(new RationalInteger(iProduct));
  }

  sort(aiPossibleNorms2.begin(), aiPossibleNorms2.end(),
       isLessThanPtrAlgebraicInteger);

  aiPossibleNorms2 = vector<AlgebraicInteger *>(
      aiPossibleNorms2.begin(),
      unique(aiPossibleNorms2.begin(), aiPossibleNorms2.end(),
             isEqualToPtrAlgebraicInteger));

  vf = new RationalInteger_VFs(aiPossibleNorms2);
}

void RationalInteger_AlVin::findVector(AlgebraicInteger *aiX0,
                                       AlgebraicInteger *aiNorm2) {
  iBilinearProducts = vector<int>(vectorsCountSecond, 0);

  findVector(dynamic_cast<RationalInteger *>(aiX0)->iVal,
             dynamic_cast<RationalInteger *>(aiNorm2)->iVal);
}

bool RationalInteger_AlVin::bCandidatePreserveEvenLattice(
    const unsigned int &iNorm2) const {
  cout << "bCandidatePreserveEvenLattice: " << iNorm2 << endl;

  long int iSumCoordinates(iVectorCurrent[0]);
  unsigned int i;
  for (i = 1; i <= dimension; i++)
    iSumCoordinates += iVectorCurrent[i];

  if (iSumCoordinates % 4 == 0)
    return true;

  iSumCoordinates *= 2;

  if (iSumCoordinates *
      (riQF[0] * iVectorCurrent[0] + riQF[1] * iVectorCurrent[1]) % 2)
    return false;

  for (i = 1; i <= dimension; i++) {
    if (iSumCoordinates *
        (-riQF[0] * iVectorCurrent[0] + riQF[i] * iVectorCurrent[i]) % 2)
      return false;
  }

  cout << "bCandidatePreserveEvenLattice: ok" << endl;

  return true;
}

void RationalInteger_AlVin::addCandidate() {
  vector<AlgebraicInteger *> aiV;
  for (auto i : iVectorCurrent)
    aiV.push_back(new RationalInteger(i));

  candidateVectors.push_back(aiV);
}

void RationalInteger_AlVin::addVectorChild(
    const std::vector<AlgebraicInteger *> &aiVector) {
  vector<int> iV;
  for (auto ai : aiVector)
    iV.push_back(dynamic_cast<RationalInteger *>(ai)->iVal);

  iVectors.push_back(iV);
}

int RationalInteger_AlVin::addVector_findWeight(
    AlgebraicInteger *aiNumerator, AlgebraicInteger *aiDenominator) {
  return -2;
}

void RationalInteger_AlVin::findVector(const unsigned int &i0,
                                       const unsigned int &iNorm2,
                                       unsigned int iIndex,
                                       unsigned int iSumComp,
                                       unsigned int iGCDComponents) {
  if (iIndex == 1) {
    iSumComp = iNorm2 + riQF[0] * i0 * i0;
    iVectorCurrent[0] = i0;
    iGCDComponents = i0;

    for (unsigned int i(0); i < vectorsCountSecond; i++)
      iBilinearProducts[i] = -riQF[0] * i0 * iVectors[i + dimension][0];
  } else if (iGCDs[iIndex] &&
             iSumComp % iGCDs[iIndex]) // equation a * x + b * y = z has a
                                       // solution only if gcd(a, b) | z
    return;

  if (!iSumComp) {
    for (; iIndex <= dimension; iIndex++)
      iVectorCurrent[iIndex] = 0;

    if (iGCDComponents == 1)
      addCandidate();

    return;
  }

  if (iIndex == dimension - 1) {
    if (bLastTwoCoefficientsQFAre1Mod4 &&
        iSumComp % 4 == 3) // Sum of two squares cannot be 3 mod 4
      return;

    unsigned int a(riQF[iIndex]), b(riQF[dimension]), c(iSumComp);

    if (iGCDs[iIndex]) {
      a /= iGCDs[iIndex];
      b /= iGCDs[iIndex];
      c *= iGCDs[iIndex];
    }

    unsigned int iGCD(ugcd(a, c));
    a /= iGCD;
    b *= iGCD;
    c /= iGCD;

    iGCD = ugcd(b, c);
    a *= iGCD;
    b /= iGCD;
    c /= iGCD;

    if (!(c % 3) && (c % 9) && jacobiSymbol(-a * b, 3) == -1)
      return;

    if (!(c % 5) && (c % 25) && jacobiSymbol(-a * b, 5) == -1)
      return;

    if (!(c % 7) && (c % 49) && jacobiSymbol(-a * b, 7) == -1)
      return;

    if (!(c % 11) && (c % 121) && jacobiSymbol(-a * b, 11) == -1)
      return;

    if (!(c % 13) && (c % 169) && jacobiSymbol(-a * b, 13) == -1)
      return;
  }

  unsigned int iMax(
      componentLessThan[iIndex]
          ? min((unsigned int)iVectorCurrent[componentLessThan[iIndex]],
                integerSqrt(iSumComp) / integerSqrt(riQF[iIndex]))
          : integerSqrt(iSumComp) / integerSqrt(riQF[iIndex]));

  iVectorCurrent[iIndex] = 0;
  if (iIndex < dimension)
    findVector(i0, iNorm2, iIndex + 1, iSumComp, iGCDComponents);

  unsigned int iLastCoeff(0);
  for (unsigned int i(ceilQuotient(iNorm2, 2 * riQF[iIndex])); i <= iMax; i++) {
    if (!(2 * i * riQF[iIndex] % iNorm2)) {
      for (unsigned int j(0); j < vectorsCountSecond; j++) {
        iBilinearProducts[j] +=
            riQF[iIndex] * (i - iLastCoeff) * iVectors[j + dimension][iIndex];

        if (iBilinearProducts[j] > 0) // iVectorCurrent is not admissible
        {
          // we restore the value of iBilinearProducts
          for (unsigned int k(0); k < vectorsCountSecond; k++)
            iBilinearProducts[k] -= riQF[iIndex] *
                                    iVectors[k + dimension][iIndex] *
                                    (k <= j ? i : iLastCoeff);

          return;
        }
      }

      iVectorCurrent[iIndex] = i;
      iLastCoeff = i;

      if (iIndex == dimension && riQF[iIndex] * i * i == iSumComp) {
        if (ugcd(i, iGCDComponents) == 1)
          addCandidate();

        for (unsigned int i(0); i < vectorsCountSecond; i++)
          iBilinearProducts[i] -=
              riQF[iIndex] * iLastCoeff * iVectors[i + dimension][iIndex];

        return;
      }

      if (iIndex < dimension) {
        if (riQF[iIndex] * i * i <= iSumComp)
          findVector(i0, iNorm2, iIndex + 1, iSumComp - riQF[iIndex] * i * i,
                     ugcd(iGCDComponents, i));
      }
    }
  }

  for (unsigned int i(0); i < vectorsCountSecond; i++)
    iBilinearProducts[i] -=
        riQF[iIndex] * iLastCoeff * iVectors[i + dimension][iIndex];
}

string RationalInteger_AlVin::get_strField() const { return string("Q"); }
