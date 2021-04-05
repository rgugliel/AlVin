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
  vectorCurrent = vector<int>(dimension + 1, 0);

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
  vector<AlgebraicInteger *> possibleNorms2(
      {new RationalInteger(1)}); ///< Possible values for (e,e)

  vector<unsigned int> primeNumbers;

  // ----------------------------------------------------------
  // Compute for ever coefficient a of the quadratic form the prime numbers
  // which divide a
  for (unsigned int i(0); i <= dimension; i++) {
    unsigned int iDivisor(3), n(riQF[i]);
    if (n == 1)
      continue;

    if (!(n % 2)) {
      primeNumbers.push_back(2);
      do {
        n /= 2;
      } while (!(n % 2));
    }

    while (n > 1) {
      if (!(n % iDivisor)) {
        primeNumbers.push_back(iDivisor);
        do {
          n /= iDivisor;
        } while (!(n % iDivisor));
      }

      iDivisor += 2;
    }
  }

  sort(primeNumbers.begin(), primeNumbers.end());
  primeNumbers = vector<unsigned int>(
      primeNumbers.begin(), unique(primeNumbers.begin(), primeNumbers.end()));

  primeNumbers.push_back(2);

  unsigned int primeNumbersCount(primeNumbers.size());
  if (primeNumbersCount > 8)
    throw(string("Number of values for (e,e) is too big"));

  // ----------------------------------------------------------
  // Powers of two
  vector<unsigned int> pow2({1});
  for (unsigned int i(1); i < primeNumbersCount; i++)
    pow2[i] = pow2[i - 1] * 2;

  // ----------------------------------------------------------
  // Compute the products of prime numbers
  unsigned int iMax(pow(2, primeNumbers.size()) - 1), j, iProduct;

  for (unsigned int i = 1; i <= iMax; i++) {
    iProduct = 1;
    for (j = 0; j < primeNumbersCount; j++) {
      if (i & pow2[j])
        iProduct *= primeNumbers[j];
    }

    possibleNorms2.push_back(new RationalInteger(iProduct));
  }

  sort(possibleNorms2.begin(), possibleNorms2.end(),
       isLessThanPtrAlgebraicInteger);

  possibleNorms2 = vector<AlgebraicInteger *>(
      possibleNorms2.begin(),
      unique(possibleNorms2.begin(), possibleNorms2.end(),
             isEqualToPtrAlgebraicInteger));

  vf = new RationalInteger_VFs(possibleNorms2);
}

void RationalInteger_AlVin::findVector(AlgebraicInteger *x0,
                                       AlgebraicInteger *norm2) {
  bilinearProducts = vector<int>(vectorsCountSecond, 0);

  findVector(dynamic_cast<RationalInteger *>(x0)->val,
             dynamic_cast<RationalInteger *>(norm2)->val);
}

bool RationalInteger_AlVin::bCandidatePreserveEvenLattice(
    const unsigned int &norm2) const {
  cout << "bCandidatePreserveEvenLattice: " << norm2 << endl;

  long int iSumCoordinates(vectorCurrent[0]);
  unsigned int i;
  for (i = 1; i <= dimension; i++)
    iSumCoordinates += vectorCurrent[i];

  if (iSumCoordinates % 4 == 0)
    return true;

  iSumCoordinates *= 2;

  if (iSumCoordinates *
      (riQF[0] * vectorCurrent[0] + riQF[1] * vectorCurrent[1]) % 2)
    return false;

  for (i = 1; i <= dimension; i++) {
    if (iSumCoordinates *
        (-riQF[0] * vectorCurrent[0] + riQF[i] * vectorCurrent[i]) % 2)
      return false;
  }

  cout << "bCandidatePreserveEvenLattice: ok" << endl;

  return true;
}

void RationalInteger_AlVin::addCandidate() {
  vector<AlgebraicInteger *> v;
  for (const auto &i : vectorCurrent)
    v.push_back(new RationalInteger(i));

  candidateVectors.push_back(v);
}

void RationalInteger_AlVin::addVectorChild(
    const std::vector<AlgebraicInteger *> &v) {
  vector<int> iV;
  for (const auto &ai : v)
    iV.push_back(dynamic_cast<RationalInteger *>(ai)->val);

  vectors.push_back(iV);
}

int RationalInteger_AlVin::addVector_findWeight(AlgebraicInteger *numerator,
                                                AlgebraicInteger *denominator) {
  return -2;
}

void RationalInteger_AlVin::findVector(const unsigned int &i0,
                                       const unsigned int &norm2,
                                       unsigned int index, unsigned int sumComp,
                                       unsigned int gcdComponents) {
  if (index == 1) {
    sumComp = norm2 + riQF[0] * i0 * i0;
    vectorCurrent[0] = i0;
    gcdComponents = i0;

    for (unsigned int i(0); i < vectorsCountSecond; i++)
      bilinearProducts[i] = -riQF[0] * i0 * vectors[i + dimension][0];
  } else if (iGCDs[index] &&
             sumComp % iGCDs[index]) // equation a * x + b * y = z has a
                                     // solution only if gcd(a, b) | z
    return;

  if (!sumComp) {
    for (; index <= dimension; index++)
      vectorCurrent[index] = 0;

    if (gcdComponents == 1)
      addCandidate();

    return;
  }

  if (index == dimension - 1) {
    if (bLastTwoCoefficientsQFAre1Mod4 &&
        sumComp % 4 == 3) // Sum of two squares cannot be 3 mod 4
      return;

    unsigned int a(riQF[index]), b(riQF[dimension]), c(sumComp);

    if (iGCDs[index]) {
      a /= iGCDs[index];
      b /= iGCDs[index];
      c *= iGCDs[index];
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
      componentLessThan[index]
          ? min((unsigned int)vectorCurrent[componentLessThan[index]],
                integerSqrt(sumComp) / integerSqrt(riQF[index]))
          : integerSqrt(sumComp) / integerSqrt(riQF[index]));

  vectorCurrent[index] = 0;
  if (index < dimension)
    findVector(i0, norm2, index + 1, sumComp, gcdComponents);

  unsigned int iLastCoeff(0);
  for (unsigned int i(ceilQuotient(norm2, 2 * riQF[index])); i <= iMax; i++) {
    if (!(2 * i * riQF[index] % norm2)) {
      for (unsigned int j(0); j < vectorsCountSecond; j++) {
        bilinearProducts[j] +=
            riQF[index] * (i - iLastCoeff) * vectors[j + dimension][index];

        if (bilinearProducts[j] > 0) // iVectorCurrent is not admissible
        {
          // we restore the value of iBilinearProducts
          for (unsigned int k(0); k < vectorsCountSecond; k++)
            bilinearProducts[k] -= riQF[index] * vectors[k + dimension][index] *
                                   (k <= j ? i : iLastCoeff);

          return;
        }
      }

      vectorCurrent[index] = i;
      iLastCoeff = i;

      if (index == dimension && riQF[index] * i * i == sumComp) {
        if (ugcd(i, gcdComponents) == 1)
          addCandidate();

        for (unsigned int i(0); i < vectorsCountSecond; i++)
          bilinearProducts[i] -=
              riQF[index] * iLastCoeff * vectors[i + dimension][index];

        return;
      }

      if (index < dimension) {
        if (riQF[index] * i * i <= sumComp)
          findVector(i0, norm2, index + 1, sumComp - riQF[index] * i * i,
                     ugcd(gcdComponents, i));
      }
    }
  }

  for (unsigned int i(0); i < vectorsCountSecond; i++)
    bilinearProducts[i] -=
        riQF[index] * iLastCoeff * vectors[i + dimension][index];
}

string RationalInteger_AlVin::get_strField() const { return string("Q"); }
