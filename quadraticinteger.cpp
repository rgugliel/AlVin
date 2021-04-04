#include "quadraticinteger.h"

int QuadraticInteger::d = 0;
int QuadraticInteger::iDiscriminant = 0;
bool QuadraticInteger::bIsOneMod4 = false;

map<unsigned int, vector<long int>> QuadraticInteger::iFundamentalUnits = {
    {2, vector<long int>({1, 1})},
    {3, vector<long int>({2, 1})},
    {5, vector<long int>({0, 1})},
    {6, vector<long int>({5, 2})},
    {7, vector<long int>({8, 3})},
    {11, vector<long int>({10, 3})},
    {13, vector<long int>({1, 1})},
    {14, vector<long int>({15, 4})},
    {17, vector<long int>({3, 2})},
    {19, vector<long int>({170, 39})},
    {21, vector<long int>({2, 1})},
    {22, vector<long int>({197, 42})},
    {23, vector<long int>({24, 5})},
    {29, vector<long int>({2, 1})},
    {31, vector<long int>({1520, 273})},
    {33, vector<long int>({19, 8})},
    {37, vector<long int>({5, 2})},
    {38, vector<long int>({37, 6})},
    {41, vector<long int>({27, 10})},
    {43, vector<long int>({3482, 531})},
    {46, vector<long int>({24335, 3588})},
    {47, vector<long int>({48, 7})},
    {51, vector<long int>({50, 7})},
    {53, vector<long int>({3, 1})},
    {57, vector<long int>({131, 40})},
    {58, vector<long int>({99, 13})},
    {59, vector<long int>({530, 69})},
    {61, vector<long int>({17, 5})},
    {62, vector<long int>({63, 8})},
    {67, vector<long int>({48842, 5967})},
    {69, vector<long int>({11, 3})},
    {71, vector<long int>({3480, 413})},
    {73, vector<long int>({943, 250})},
    {74, vector<long int>({43, 5})},
    {77, vector<long int>({4, 1})},
    {78, vector<long int>({53, 6})},
    {79, vector<long int>({80, 9})},
    {82, vector<long int>({9, 1})},
    {83, vector<long int>({82, 9})},
    {85, vector<long int>({4, 1})},
    {86, vector<long int>({10405, 1122})},
    {89, vector<long int>({447, 106})},
    {91, vector<long int>({1574, 165})},
    {93, vector<long int>({13, 3})},
    {94, vector<long int>({2143295, 221064})},
    {95, vector<long int>({39, 4})},
    {97, vector<long int>({5035, 1138})}};

map<unsigned int, vector<unsigned int>> QuadraticInteger::iPellMinimalSolution =
    {{2, vector<unsigned int>({3, 2})},
     {3, vector<unsigned int>({2, 1})},
     {5, vector<unsigned int>({9, 4})},
     {6, vector<unsigned int>({5, 2})},
     {7, vector<unsigned int>({8, 3})},
     {11, vector<unsigned int>({10, 3})},
     {13, vector<unsigned int>({649, 180})},
     {14, vector<unsigned int>({15, 4})},
     {17, vector<unsigned int>({33, 8})},
     {19, vector<unsigned int>({170, 39})},
     {21, vector<unsigned int>({55, 12})},
     {22, vector<unsigned int>({197, 42})},
     {23, vector<unsigned int>({24, 5})},
     {29, vector<unsigned int>({9801, 1820})},
     {31, vector<unsigned int>({1520, 273})},
     {33, vector<unsigned int>({23, 4})},
     {37, vector<unsigned int>({73, 12})},
     {38, vector<unsigned int>({37, 6})},
     {41, vector<unsigned int>({2049, 320})},
     {43, vector<unsigned int>({3482, 531})},
     {46, vector<unsigned int>({24335, 3588})},
     {47, vector<unsigned int>({48, 7})},
     {51, vector<unsigned int>({50, 7})},
     {53, vector<unsigned int>({66249, 9100})},
     {57, vector<unsigned int>({151, 20})},
     {58, vector<unsigned int>({19603, 2574})},
     {59, vector<unsigned int>({530, 69})},
     {61, vector<unsigned int>({1766319049, 226153980})},
     {62, vector<unsigned int>({63, 8})},
     {67, vector<unsigned int>({48842, 5967})},
     {69, vector<unsigned int>({7775, 936})},
     {71, vector<unsigned int>({3480, 413})},
     {73, vector<unsigned int>({2281249, 267000})},
     {74, vector<unsigned int>({3699, 430})},
     {77, vector<unsigned int>({351, 40})},
     {78, vector<unsigned int>({53, 6})},
     {79, vector<unsigned int>({80, 9})},
     {82, vector<unsigned int>({163, 18})},
     {83, vector<unsigned int>({82, 9})},
     {85, vector<unsigned int>({285769, 30996})},
     {86, vector<unsigned int>({10405, 1122})},
     {89, vector<unsigned int>({500001, 53000})},
     {91, vector<unsigned int>({1574, 165})},
     {93, vector<unsigned int>({12151, 1260})},
     {94, vector<unsigned int>({2143295, 221064})},
     {95, vector<unsigned int>({39, 4})},
     {97, vector<unsigned int>({62809633, 6377352})}};

QuadraticInteger::QuadraticInteger() : a(0), b(0) {}

QuadraticInteger::QuadraticInteger(const int &iVal) : a(iVal), b(0) {}

QuadraticInteger::QuadraticInteger(const QuadraticInteger &qi)
    : a(qi.a), b(qi.b) {}

QuadraticInteger::QuadraticInteger(int a, int b) : a(a), b(b) {}
bool QuadraticInteger::bIsDAdmissible(const unsigned int &d) {
  return (iPellMinimalSolution.find(d) != iPellMinimalSolution.end());
}

void QuadraticInteger::set_d(const unsigned int &dnew) {
  d = (int)dnew;
  d = iRemoveSquareFactors(d);

  if (iPellMinimalSolution.find(d) == iPellMinimalSolution.end())
    throw(string("Invalid value for d: " + std::to_string(d)));

  bIsOneMod4 = (d % 4 == 1);
  iDiscriminant = bIsOneMod4 ? d : 4 * d;
}

QuadraticInteger::~QuadraticInteger() {}

int QuadraticInteger::iValuation(const QuadraticInteger &qi) {
  if (a == 0 && b == 0)
    return -1;

  int iVal(0);
  QuadraticInteger qiTemp(*this);

  while (qiTemp.divideByIfDivisible(&qi))
    iVal++;

  return iVal;
}

vector<QuadraticInteger> QuadraticInteger::qiPrimeFactors() const {
  if ((a == 0 && b == 0) || isInvertible())
    return vector<QuadraticInteger>(0);

  vector<QuadraticInteger> qiPrimesFactors;
  auto iNormPrimesFactors(primeFactors((unsigned long)abs(b ? iNorm() : a)));

  for (auto iNPF : iNormPrimesFactors) {
    auto qiTemp(QuadraticInteger::qiFactorsRationalPrime(iNPF));

    for (auto qiF : qiTemp) {
      if (isDivisibleBy(
              &qiF)) // qiF divides either *this or the conjugate of *this
        qiPrimesFactors.push_back(qiF);
    }
  }

  return qiPrimesFactors;
}

std::map<QuadraticInteger, unsigned int>
QuadraticInteger::qiPrimeDecomposition() const {
  // TODO: mettre en cache?

  if ((a == 0 && b == 0) || isInvertible())
    return map<QuadraticInteger, unsigned int>();

  map<QuadraticInteger, unsigned int> qiDecomp;
  vector<unsigned long> iDecom(
      primeFactors<unsigned long>(abs(b ? iNorm() : a)));
  unsigned int iPower;

  QuadraticInteger qiTemp(*this);

  for (auto it : iDecom) {
    vector<QuadraticInteger> qiFactors(qiFactorsRationalPrime(it));
    for (auto itQI : qiFactors) {
      iPower = 0;
      while (qiTemp.divideByIfDivisible(&itQI))
        iPower++;

      if (iPower)
        qiDecomp[itQI] = iPower;
    }
  }

  return qiDecomp;
}

vector<QuadraticInteger>
QuadraticInteger::qiFactorsRationalPrime(const unsigned int &iPrime,
                                         bool bWithMultiplicities) {
  // TODO: mettre en cache?
  if (!d)
    throw(string("d must be specified"));

  if (iPrime % 2 && d % iPrime) {
    if (jacobiSymbol(d, iPrime) == -1)
      return vector<QuadraticInteger>({QuadraticInteger(iPrime)});

    if (bIsOneMod4) {
      array<long int, 2> iCoeffs(iPellEquation(4 * iPrime));

      vector<QuadraticInteger> qiFactors(
          2, QuadraticInteger((iCoeffs[0] - iCoeffs[1]) / 2, iCoeffs[1]));
      qiFactors[1].conjugate();

      return qiFactors;
    } else {
      array<long int, 2> iCoeffs(iPellEquation(iPrime));

      vector<QuadraticInteger> qiFactors(
          2, QuadraticInteger(iCoeffs[0], iCoeffs[1]));
      qiFactors[1].conjugate();

      return qiFactors;
    }
  } else if (d % iPrime == 0) // p divides d
  {
    if (bIsOneMod4) {
      array<long int, 2> iCoeffs(iPellEquation(4 * iPrime));

      return vector<QuadraticInteger>(
          (int)bWithMultiplicities + 1,
          QuadraticInteger((iCoeffs[0] - iCoeffs[1]) / 2, iCoeffs[1]));
    } else {
      array<long int, 2> iCoeffs(iPellEquation(iPrime));

      return vector<QuadraticInteger>((int)bWithMultiplicities + 1,
                                      QuadraticInteger(iCoeffs[0], iCoeffs[1]));
    }
  } else // p=2 and 2 does not divide d
  {
    if (d % 4 == 3) {
      array<long int, 2> iCoeffs(iPellEquation(iPrime));

      return vector<QuadraticInteger>((int)bWithMultiplicities + 1,
                                      QuadraticInteger(iCoeffs[0], iCoeffs[1]));
    } else if (d % 8 == 1) {
      array<long int, 2> iCoeffs(iPellEquation(4 * iPrime));

      vector<QuadraticInteger> qiFactors(
          2, QuadraticInteger((iCoeffs[0] - iCoeffs[1]) / 2, iCoeffs[1]));
      qiFactors[1].conjugate();

      return qiFactors;
    } else // d % 8 == 5
      return vector<QuadraticInteger>({QuadraticInteger(iPrime)});
  }

  throw(string("qiFactorsRationalPrime: unknown case (d=" + std::to_string(d) +
               ", p=" + std::to_string(iPrime) + ")"));
}

// TODO: optimiser?
array<long int, 2> QuadraticInteger::iPellEquation(const unsigned int &iN) {
  long int iXTemp, iSqrt;
  if (iN == (unsigned int)d)
    return array<long int, 2>({0, 1});

  for (unsigned int iNegative(0); iNegative <= 1; iNegative++) {
    unsigned int L1, L2;

    // TODO: regarder si iSQRTsup
    if (iNegative == 0) {
      L1 = 0;
      L2 = integerSqrt(((iN * (iPellMinimalSolution[d][0] - 1)) % (2 * d))
                           ? (iN * (iPellMinimalSolution[d][0] - 1)) / (2 * d) +
                                 1
                           : (iN * (iPellMinimalSolution[d][0] - 1)) / (2 * d));
    } else {
      L1 = integerSqrt(iN / d);
      L2 = integerSqrt(((iN * (iPellMinimalSolution[d][0] + 1)) % (2 * d))
                           ? (iN * (iPellMinimalSolution[d][0] + 1)) / (2 * d) +
                                 1
                           : (iN * (iPellMinimalSolution[d][0] + 1)) / (2 * d));
    }

    for (unsigned int y(L1); y <= L2; y++) {
      iXTemp = (iNegative ? -1 : 1) * iN + d * y * y;

      if (iXTemp >= 0) // if it is a square
      {
        iSqrt = integerSqrt((unsigned long int)iXTemp);

        if (iSqrt * iSqrt == iXTemp)
          return array<long int, 2>({iSqrt, y});
      }
    }
  }

  throw(string("Pell equation was not solved for d=" + std::to_string(d) +
               " and n=" + std::to_string(iN)));
}

double QuadraticInteger::to_double() const {
  if (bIsOneMod4)
    return (a + (1 + sqrt(d)) * b / 2);
  else
    return (a + b * sqrt(d));
}

bool QuadraticInteger::isEqualTo(const int &n) const {
  return (a == n && b == 0);
}

bool QuadraticInteger::isEqualTo(const AlgebraicInteger &ai) const {
  QuadraticInteger qi(dynamic_cast<const QuadraticInteger &>(ai));

  return (a == qi.a && b == qi.b);
}

bool QuadraticInteger::isGreaterThan(const int &n) const {
  long int aTemp(bIsOneMod4 ? (2 * (n - a) - b) : (n - a));
  long int bTemp(-b);

  if (aTemp < 0) {
    if (bTemp <= 0)
      return true;
    else
      return (bTemp * bTemp * d < aTemp * aTemp);
  } else {
    if (bTemp >= 0)
      return false;
    else
      return (aTemp * aTemp < bTemp * bTemp * d);
  }
}

bool QuadraticInteger::bIsGreaterThan(const long int &n) const {
  long int aTemp(bIsOneMod4 ? (2 * (n - a) - b) : (n - a));
  long int bTemp(-b);

  if (aTemp < 0) {
    if (bTemp <= 0)
      return true;
    else
      return (bTemp * bTemp * d < aTemp * aTemp);
  } else {
    if (bTemp >= 0)
      return false;
    else
      return (aTemp * aTemp < bTemp * bTemp * d);
  }
}

bool QuadraticInteger::isGreaterOEThan(const int &n) const {
  return ((b == 0 && a == n) || isGreaterThan(n));
}

bool operator<(const QuadraticInteger &qi1, const QuadraticInteger &qi2) {
  return qi1.isLessThan(qi2);
}

bool QuadraticInteger::isLessThan(const AlgebraicInteger &ai) const {
  long int x(dynamic_cast<const QuadraticInteger &>(ai).a);
  long int y(dynamic_cast<const QuadraticInteger &>(ai).b);

  long int aTemp(bIsOneMod4 ? (2 * (a - x) + b - y) : a - x);
  long int bTemp(b - y);

  if (aTemp < 0) {
    if (bTemp <= 0)
      return true;
    else
      return (bTemp * bTemp * d < aTemp * aTemp);
  } else {
    if (bTemp >= 0)
      return false;
    else
      return (aTemp * aTemp < bTemp * bTemp * d);
  }
}

bool QuadraticInteger::isLessOEThan(const AlgebraicInteger &ai) const {
  long int x(dynamic_cast<const QuadraticInteger &>(ai).a);
  long int y(dynamic_cast<const QuadraticInteger &>(ai).b);

  if (x == a && y == b)
    return true;

  long int aTemp(bIsOneMod4 ? (2 * (a - x) + b - y) : a - x);
  long int bTemp(b - y);

  if (aTemp < 0) {
    if (bTemp <= 0)
      return true;
    else
      return (bTemp * bTemp * d < aTemp * aTemp);
  } else {
    if (bTemp >= 0)
      return false;
    else
      return (aTemp * aTemp < bTemp * bTemp * d);
  }
}

bool QuadraticInteger::isLessThan(const int &n) const {
  long int aTemp(bIsOneMod4 ? (2 * (a - n) + b) : a - n);

  if (aTemp < 0) {
    if (b <= 0)
      return true;
    else
      return (b * b * d < aTemp * aTemp);
  } else {
    if (b >= 0)
      return false;
    else
      return (aTemp * aTemp < b * b * d);
  }
}

bool QuadraticInteger::bIsLessThan(const long int &n) const {
  long int aTemp(bIsOneMod4 ? (2 * (a - n) + b) : a - n);

  if (aTemp < 0) {
    if (b <= 0)
      return true;
    else
      return (b * b * d < aTemp * aTemp);
  } else {
    if (b >= 0)
      return false;
    else
      return (aTemp * aTemp < b * b * d);
  }
}

void QuadraticInteger::opp() {
  a *= -1;
  b *= -1;
}

void QuadraticInteger::conjugate() {
  if (bIsOneMod4)
    a += b;

  b *= -1;
}

void QuadraticInteger::substract(const AlgebraicInteger *ai) {
  a -= dynamic_cast<const QuadraticInteger *>(ai)->a;
  b -= dynamic_cast<const QuadraticInteger *>(ai)->b;
}

void QuadraticInteger::add(const AlgebraicInteger *ai) {
  a += dynamic_cast<const QuadraticInteger *>(ai)->a;
  b += dynamic_cast<const QuadraticInteger *>(ai)->b;
}

void QuadraticInteger::multiplyBy(const int &n) {
  a *= n;
  b *= n;
}

void QuadraticInteger::multiplyBy(const AlgebraicInteger *ai) {
  long int iA1(a), iB1(b);
  long int iA2(dynamic_cast<const QuadraticInteger *>(ai)->a);
  long int iB2(dynamic_cast<const QuadraticInteger *>(ai)->b);

  if (bIsOneMod4) {
    a = iA1 * iA2 + iB1 * iB2 * (d - 1) / 4;
    b = iB1 * iA2 + iA1 * iB2 + iB1 * iB2;
  } else {
    a = iA1 * iA2 + d * iB1 * iB2;
    b = iA1 * iB2 + iB1 * iA2;
  }
}

long int QuadraticInteger::floor() const {
  if (b == 0)
    return a;

  long int iRes(to_double());

  while (bIsGreaterThan(iRes + 1))
    iRes++;

  while (bIsLessThan(iRes))
    iRes--;

  return iRes;
}

bool QuadraticInteger::divideByIfDivisible(const AlgebraicInteger *ai) {
  const QuadraticInteger *qi(dynamic_cast<const QuadraticInteger *>(ai));

  long int iA1(a), iB1(b);
  long int iA2(dynamic_cast<const QuadraticInteger *>(ai)->a);
  long int iB2(dynamic_cast<const QuadraticInteger *>(ai)->b);

  long int iNorm2(qi->iNorm());

  if (this->iNorm() % iNorm2)
    return false;

  if (bIsOneMod4) {
    if ((a * (2 * iA2 + iB2) + b * iA2 - b * iB2 * (d - 1) / 2) % iNorm2)
      return false;

    a = (iA1 * (iA2 + iB2) - iB1 * iB2 * (d - 1) / 4) / iNorm2;
    b = (iB1 * iA2 - iA1 * iB2) / iNorm2;
  } else {
    if (((a * iA2 - b * iB2 * d) * 2) % iNorm2)
      return false;

    a = (iA1 * iA2 - iB1 * iB2 * d) / iNorm2;
    b = (iB1 * iA2 - iA1 * iB2) / iNorm2;
  }

  return true;
}

void QuadraticInteger::divideBy(const AlgebraicInteger *ai) {
  long int iA1(a), iB1(b);
  long int iA2(dynamic_cast<const QuadraticInteger *>(ai)->a);
  long int iB2(dynamic_cast<const QuadraticInteger *>(ai)->b);
  long int iNorm2(dynamic_cast<const QuadraticInteger *>(ai)->iNorm());

  if (bIsOneMod4) {
    a = (iA1 * (iA2 + iB2) - iB1 * iB2 * (d - 1) / 4) / iNorm2;
    b = (iB1 * iA2 - iA1 * iB2) / iNorm2;
  } else {
    a = (iA1 * iA2 - iB1 * iB2 * d) / iNorm2;
    b = (iB1 * iA2 - iA1 * iB2) / iNorm2;
  }
}

bool QuadraticInteger::isDivisibleBy(const AlgebraicInteger *ai) const {
  const QuadraticInteger *qi(dynamic_cast<const QuadraticInteger *>(ai));

  long int iNorm2(qi->iNorm());

  if (this->iNorm() % iNorm2)
    return false;

  if (bIsOneMod4) {
    if ((a * (2 * qi->a + qi->b) + b * qi->a - b * qi->b * (d - 1) / 2) %
        iNorm2)
      return false;
  } else {
    if (((a * qi->a - b * qi->b * d) * 2) % iNorm2)
      return false;
  }

  return true;
}

bool QuadraticInteger::bIsAssociateTo(QuadraticInteger qi2) {
  if (qi2.divideByIfDivisible(this)) {
    if (qi2.isInvertible())
      return true;
  }

  return false;
}

void QuadraticInteger::gcd(const AlgebraicInteger *ai) {
  const QuadraticInteger *qi(dynamic_cast<const QuadraticInteger *>(ai));

  // -----------------------------------------
  // If one of the two is zero
  if (ai->isEqualTo(0))
    return;

  if (this->isEqualTo(0)) {
    const QuadraticInteger *qi(dynamic_cast<const QuadraticInteger *>(ai));
    a = qi->a;
    b = qi->b;
    return;
  }

  // -----------------------------------------
  // If one is bIsInvertible
  unsigned long int iNorm1(abs(this->iNorm())), iNorm2(abs(qi->iNorm()));
  if (iNorm1 == 1 || iNorm2 == 1) {
    a = 1;
    b = 0;
    return;
  }

  // -----------------------------------------
  // Some other trivial cases
  if (a == qi->a && b == qi->b)
    return;
  else if (qi->isDivisibleBy(this))
    return;
  else if (this->isDivisibleBy(qi)) {
    a = qi->a;
    b = qi->b;
    return;
  }

  // -----------------------------------------
  // Computations
  unsigned long int iGCD(ugcd<unsigned long int>(iNorm1, iNorm2));
  if (iGCD == 1) {
    a = 1;
    b = 0;
    return;
  }

  vector<QuadraticInteger> qiGCD(QuadraticInteger(iGCD).qiPrimeFactors());
  QuadraticInteger qi1(*this), qi2(*qi);

  a = 1;
  b = 0;

  for (auto qiFactor : qiGCD) {
    while (qi1.divideByIfDivisible(&qiFactor) &&
           qi2.divideByIfDivisible(&qiFactor))
      this->multiplyBy(&qiFactor);
  }
}

bool QuadraticInteger::isSquareOfIvertible() const {
  if (this->isInvertible()) {
    QuadraticInteger qiUnit(QuadraticInteger::iFundamentalUnits[5][0],
                            QuadraticInteger::iFundamentalUnits[5][1]);
    QuadraticInteger qiTemp(*this);
    unsigned int iPower(0);

    while (!qiTemp.isEqualTo(1)) {
      qiTemp.divideBy(&qiUnit);
      iPower++;
    }

    return (iPower % 2 == 0);
  }
  return false;
}

bool QuadraticInteger::isInvertible() const {
  long int iNorm(this->iNorm());
  return (iNorm == 1 || iNorm == -1);
}

void QuadraticInteger::removeSquareFactors() {
  map<QuadraticInteger, unsigned int> qiDecomp(this->qiPrimeDecomposition());
  unsigned int i;
  for (auto it : qiDecomp) {
    unsigned int iMax(it.second % 2 ? it.second - 1 : it.second);
    for (i = 1; i <= iMax; i++)
      this->divideBy(&it.first);
  }
}

void QuadraticInteger::set(AlgebraicInteger *ai) {
  QuadraticInteger *qi(dynamic_cast<QuadraticInteger *>(ai));

  a = qi->a;
  b = qi->b;
}

void QuadraticInteger::set(const int &n) {
  a = n;
  b = 0;
}

long int QuadraticInteger::iNorm() const {
  return (bIsOneMod4 ? (a * (a + b) - b * b * (d - 1) / 4)
                     : (a * a - b * b * (long int)d));
}

long int QuadraticInteger::iTrace() const {
  return (bIsOneMod4 ? (2 * a + b) : 2 * a);
}

long int QuadraticInteger::iSQRTsup_quotient(const QuadraticInteger &qiNum,
                                             const QuadraticInteger &qiDen) {
  long int x, y, z, s;
  long int iNumTrace(qiNum.iTrace()), iDenTrace(qiDen.iTrace());

  if (bIsOneMod4) {
    x = iNumTrace * iDenTrace - qiNum.b * qiDen.b * d;
    y = qiNum.b * iDenTrace - qiDen.b * iNumTrace;
    z = 4 * qiDen.iNorm();
  } else {
    x = qiNum.a * qiDen.a - qiNum.b * qiDen.b * d;
    y = qiNum.b * qiDen.a - qiNum.a * qiDen.b;
    z = qiDen.iNorm();
  }

  if (z < 0) {
    y *= -1;
    z *= -1;
    x *= -1;
  }

  if (y >= 0) {
    s = sqrtSupQuotient<unsigned long int>(
        x + sqrtSup<unsigned long int>(y * y * d), z);

    while ((z * (s + 1) * (s + 1) - x) * (z * (s + 1) * (s + 1) - x) <=
           y * y * d)
      s--;
  } else {
    long int is(integerSqrt((unsigned long int)y * y * d));
    s = sqrtSupQuotient<unsigned long int>(x - is, z);

    long int iTemp(x - (s - 1) * (s - 1) * z);

    while (iTemp < 0 || iTemp * iTemp <= is) {
      s--;
      iTemp = x - (s - 1) * (s - 1) * z;
    }
  }

  return s;
}

long int QuadraticInteger::iSQRT_quotient(const QuadraticInteger &qiNum,
                                          const QuadraticInteger &qiDen) {
  long int x, y, z, s;
  long int iNumTrace(qiNum.iTrace()), iDenTrace(qiDen.iTrace());

  if (bIsOneMod4) {
    x = iNumTrace * iDenTrace - qiNum.b * qiDen.b * d;
    y = qiNum.b * iDenTrace - qiDen.b * iNumTrace;
    z = 4 * qiDen.iNorm();
  } else {
    x = qiNum.a * qiDen.a - qiNum.b * qiDen.b * d;
    y = qiNum.b * qiDen.a - qiNum.a * qiDen.b;
    z = qiDen.iNorm();
  }

  if (z < 0) {
    y *= -1;
    z *= -1;
    x *= -1;
  }

  if (y >= 0) {
    s = sqrtQuotient<unsigned long int>(
        x + integerSqrt((unsigned long int)y * y * d), z);

    while ((z * (s + 1) * (s + 1) - x) * (z * (s + 1) * (s + 1) - x) <=
           y * y * d)
      s++;
  } else {
    long int is(integerSqrt((unsigned long int)y * y * d));
    s = sqrtSupQuotient<unsigned long int>(x - is, z);

    long int iTemp(z * s * s - x);

    while ((iTemp > 0) || (iTemp <= 0 && iTemp * iTemp < y * y * d)) {
      s--;
      iTemp = z * s * s - x;
    }
  }

  return s;
}

AlgebraicInteger *QuadraticInteger::copy() const {
  return new QuadraticInteger(*this);
}

AlgebraicInteger *QuadraticInteger::copyToInteger(const int &n) const {
  return new QuadraticInteger(n);
}

std::ostream &QuadraticInteger::print(std::ostream &o) const {
  if (a == 0) {
    if (b == 0)
      o << "0";
    else
      o << b << " * T"
        << "(" << d << ")";
  } else {
    if (b == 0)
      o << a;
    else
      o << a << (b < 0 ? " - " : " + ") << abs(b) << " * T"
        << "(" << d << ")";
  }

  return o;
}

std::string QuadraticInteger::to_string(const string &strFormat,
                                        const bool &bProtect) const {
  if (strFormat == "latex") {
    if (bIsOneMod4) {
      if (a == 0) {
        if (b == 0)
          return "0";
        else
          return (abs(b) == 1 ? (b == 1 ? "" : "-")
                              : std::to_string(b) + "\\cdot") +
                 "\\frac{1+\\sqrt{" + std::to_string(d) + "}}{2}";
      } else {
        if (b == 0)
          return std::to_string(a);
        else
          return std::to_string(a) + (b < 0 ? "-" : "+") +
                 (abs(b) == 1 ? "" : std::to_string(abs(b)) + "\\cdot") +
                 "\\frac{1+\\sqrt{" + std::to_string(d) + "}}{2}";
      }
    } else {
      if (a == 0) {
        if (b == 0)
          return "0";
        else
          return (abs(b) == 1 ? (b == 1 ? "" : "-")
                              : std::to_string(b) + "\\cdot") +
                 "\\sqrt{" + std::to_string(d) + "}";
      } else {
        if (b == 0)
          return std::to_string(a);
        else
          return std::to_string(a) + (b < 0 ? "-" : "+") +
                 (abs(b) == 1 ? "" : std::to_string(abs(b)) + "\\cdot") +
                 "\\sqrt{" + std::to_string(d) + "}";
      }
    }
  } else if (strFormat == "mathematica") {
    if (bIsOneMod4) {
      if (a == 0) {
        if (b == 0)
          return "0";
        else
          return std::to_string(b) + "*(1+Sqrt[" + std::to_string(d) + "])/2";
      } else {
        if (b == 0)
          return std::to_string(a);
        else
          return (bProtect ? "(" : "") + std::to_string(a) +
                 (b < 0 ? "-" : "+") + std::to_string(abs(b)) + "*(1+Sqrt[" +
                 std::to_string(d) + (bProtect ? "])/2)" : "])/2");
      }
    } else {
      if (a == 0) {
        if (b == 0)
          return "0";
        else
          return std::to_string(b) + "*Sqrt[" + std::to_string(d) + "]";
      } else {
        if (b == 0)
          return std::to_string(a);
        else
          return (bProtect ? "(" : "") + std::to_string(a) +
                 (b < 0 ? "-" : "+") + std::to_string(abs(b)) + "*Sqrt[" +
                 std::to_string(d) + (bProtect ? "])" : "]");
      }
    }
  } else if (strFormat == "pari") {
    if (bIsOneMod4) {
      if (a == 0) {
        if (b == 0)
          return "0";
        else
          return std::to_string(b) + "*(1+sqrt(" + std::to_string(d) + "))/2";
      } else {
        if (b == 0)
          return std::to_string(a);
        else
          return (bProtect ? "(" : "") + std::to_string(a) +
                 (b < 0 ? "-" : "+") + std::to_string(abs(b)) + "*(1+sqrt(" +
                 std::to_string(d) + (bProtect ? "))/2)" : "))/2");
      }
    } else {
      if (a == 0) {
        if (b == 0)
          return "0";
        else
          return std::to_string(b) + "*sqrt(" + std::to_string(d) + ")";
      } else {
        if (b == 0)
          return std::to_string(a);
        else
          return (bProtect ? "(" : "") + std::to_string(a) +
                 (b < 0 ? "-" : "+") + std::to_string(abs(b)) + "*sqrt(" +
                 std::to_string(d) + (bProtect ? "))" : ")");
      }
    }

  } else if (strFormat == "filename") {
    if (a == 0) {
      if (b == 0)
        return "0";
      else
        return (b == 1 ? "" : std::to_string(b)) + "T" + std::to_string(d);
    } else {
      if (b == 0)
        return std::to_string(a);
      else
        return (bProtect ? "(" : "") + std::to_string(a) + (b < 0 ? "-" : "+") +
               (abs(b) == 1 ? "" : std::to_string(abs(b))) + "T" +
               std::to_string(d) + (bProtect ? ")" : "");
    }
  } else {
    if (a == 0) {
      if (b == 0)
        return "0";
      else
        return std::to_string(b) + "*T" + std::to_string(d);
    } else {
      if (b == 0)
        return std::to_string(a);
      else
        return (bProtect ? "(" : "") + std::to_string(a) + (b < 0 ? "-" : "+") +
               std::to_string(abs(b)) + "*T" + std::to_string(d) +
               (bProtect ? ")" : "");
    }
  }
}

string QuadraticInteger::get_classname() const { return "QuadraticInteger"; }

// ------------------------------------------------------------------------
// Operators
QuadraticInteger &QuadraticInteger::operator=(const QuadraticInteger &qi) {
  a = qi.a;
  b = qi.b;

  return *this;
}

QuadraticInteger &QuadraticInteger::operator/=(const QuadraticInteger &qi) {
  this->divideBy(&qi);
  return *this;
}

QuadraticInteger &QuadraticInteger::operator*=(const QuadraticInteger &qi) {
  this->multiplyBy(&qi);
  return *this;
}

QuadraticInteger QuadraticInteger::operator+(const QuadraticInteger &qi) const {
  return QuadraticInteger(a + qi.a, b + qi.b);
}

QuadraticInteger QuadraticInteger::operator*(const QuadraticInteger &qi) const {
  QuadraticInteger qiRes(*this);
  qiRes.multiplyBy(&qi);

  return qiRes;
}

QuadraticInteger QuadraticInteger::operator-(const QuadraticInteger &qi) const {
  return QuadraticInteger(a - qi.a, b - qi.b);
}

QuadraticInteger QuadraticInteger::operator-() const {
  return QuadraticInteger(-a, -b);
}

bool QuadraticInteger::operator>(const QuadraticInteger &ri) {
  return !isLessOEThan(ri);
}
