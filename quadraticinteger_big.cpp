#include "quadraticinteger_big.h"

int QuadraticIntegerBig::d = 0;
mpf_class QuadraticIntegerBig::sqrtd = 0;
int QuadraticIntegerBig::iDiscriminant = 0;
bool QuadraticIntegerBig::bIsOneMod4 = false;

map<unsigned int, vector<long int>> QuadraticIntegerBig::iFundamentalUnits = {
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

map<unsigned int, vector<unsigned int>>
    QuadraticIntegerBig::iPellMinimalSolution = {
        {2, vector<unsigned int>({3, 2})},
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

QuadraticIntegerBig::QuadraticIntegerBig() : a(0), b(0) {}

QuadraticIntegerBig::QuadraticIntegerBig(const int &iVal) : a(iVal), b(0) {}

QuadraticIntegerBig::QuadraticIntegerBig(const mpz_class &iVal)
    : a(iVal), b(0) {}

QuadraticIntegerBig::QuadraticIntegerBig(const QuadraticIntegerBig &qi)
    : a(qi.a), b(qi.b) {}

QuadraticIntegerBig::QuadraticIntegerBig(const QuadraticInteger &qi)
    : a(qi.a), b(qi.b) {}

QuadraticIntegerBig::QuadraticIntegerBig(const int &a, const int &b)
    : a(a), b(b) {}

QuadraticIntegerBig::QuadraticIntegerBig(const mpz_class &a, const mpz_class &b)
    : a(a), b(b) {}

bool QuadraticIntegerBig::bIsDAdmissible(const unsigned int &d) {
  return (iPellMinimalSolution.find(d) != iPellMinimalSolution.end());
}

void QuadraticIntegerBig::set_d(const unsigned int &dnew) {
  d = (int)dnew;
  d = iRemoveSquareFactors(d);

  if (iPellMinimalSolution.find(d) == iPellMinimalSolution.end())
    throw(string("Invalid value for d: " + std::to_string(d)));

  bIsOneMod4 = (d % 4 == 1);
  iDiscriminant = bIsOneMod4 ? d : 4 * d;
  sqrtd = sqrt(d);
}

QuadraticIntegerBig::~QuadraticIntegerBig() {}

int QuadraticIntegerBig::iValuation(const QuadraticIntegerBig &qi) {
  if (a == 0 && b == 0)
    return -1;

  int iVal(0);
  QuadraticIntegerBig qiTemp(*this);

  while (qiTemp.divideByIfDivisible(&qi))
    iVal++;

  return iVal;
}

vector<QuadraticIntegerBig> QuadraticIntegerBig::qiPrimeFactors() const {
  if ((a == 0 && b == 0) || isInvertible())
    return vector<QuadraticIntegerBig>(0);

  vector<QuadraticIntegerBig> qiPrimesFactors;

  mpz_class iN(abs(iNorm()));
  if (!iN.fits_ulong_p()) {
    cout << *this << endl;
    cout << iN << endl;
    throw(string("QuadraticIntegerBig::qiPrimeFactors(): Norm too big"));
  }
  auto iNormPrimesFactors(primeFactors(iN.get_ui()));

  for (auto iNPF : iNormPrimesFactors) {
    auto qiTemp(QuadraticIntegerBig::qiFactorsRationalPrime(iNPF));

    for (auto qiF : qiTemp) {
      if (isDivisibleBy(
              &qiF)) // qiF divides either *this or the conjugate of *this
        qiPrimesFactors.push_back(qiF);
    }
  }

  return qiPrimesFactors;
}

array<long int, 2> QuadraticIntegerBig::iPellEquation(const unsigned int &iN) {
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

        if (iSqrt * iSqrt == iXTemp &&
            (!bIsOneMod4 ||
             ((iXTemp % 2) == (y % 2)))) // if iXTemp is a square and if the x,y
                                         // have the same parity if d=1(4)
          return array<long int, 2>({iSqrt, y});
      }
    }
  }

  throw(string("Pell equation was not solved for d=" + std::to_string(d) +
               " and n=" + std::to_string(iN)));
}

std::map<QuadraticIntegerBig, unsigned int>
QuadraticIntegerBig::qiPrimeDecomposition() const {
  // TODO: mettre en cache?

  if ((a == 0 && b == 0) || isInvertible())
    return map<QuadraticIntegerBig, unsigned int>();

  mpz_class iN(abs(iNorm()));
  if (!iN.fits_ulong_p())
    throw(string("QuadraticIntegerBig::qiPrimeDecomposition(): Norm too big"));

  map<QuadraticIntegerBig, unsigned int> qiDecomp;
  map<unsigned long, unsigned int> iDecom(
      primeDecomposition<long unsigned>(iN.get_ui()));
  unsigned int iPower;

  QuadraticIntegerBig qiTemp(*this);

  for (auto it : iDecom) {
    vector<QuadraticIntegerBig> qiFactors(qiFactorsRationalPrime(it.first));
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

vector<QuadraticIntegerBig>
QuadraticIntegerBig::qiFactorsRationalPrime(const unsigned int &iPrime,
                                            bool bWithMultiplicities) {
  // TODO: mettre en cache?
  if (!d)
    throw(string("d must be specified"));

  if (iPrime % 2 && d % iPrime) {
    if (jacobiSymbol(d, iPrime) == -1)
      return vector<QuadraticIntegerBig>({QuadraticIntegerBig(iPrime)});

    if (bIsOneMod4) {
      array<long int, 2> iCoeffs(
          QuadraticIntegerBig::iPellEquation(4 * iPrime));

      vector<QuadraticIntegerBig> qiFactors(
          2, QuadraticIntegerBig((iCoeffs[0] - iCoeffs[1]) / 2, iCoeffs[1]));
      qiFactors[1].conjugate();

      return qiFactors;
    } else {
      array<long int, 2> iCoeffs(QuadraticIntegerBig::iPellEquation(iPrime));

      vector<QuadraticIntegerBig> qiFactors(
          2, QuadraticIntegerBig(iCoeffs[0], iCoeffs[1]));
      qiFactors[1].conjugate();

      return qiFactors;
    }
  } else if (d % iPrime == 0) // p divides d
  {
    if (bIsOneMod4) {
      array<long int, 2> iCoeffs(
          QuadraticIntegerBig::iPellEquation(4 * iPrime));

      return vector<QuadraticIntegerBig>(
          (int)bWithMultiplicities + 1,
          QuadraticIntegerBig((iCoeffs[0] - iCoeffs[1]) / 2, iCoeffs[1]));
    } else {
      array<long int, 2> iCoeffs(QuadraticIntegerBig::iPellEquation(iPrime));

      return vector<QuadraticIntegerBig>(
          (int)bWithMultiplicities + 1,
          QuadraticIntegerBig(iCoeffs[0], iCoeffs[1]));
    }
  } else // p=2 and 2 does not divide d
  {
    if (d % 4 == 3) {
      array<long int, 2> iCoeffs(QuadraticIntegerBig::iPellEquation(iPrime));

      return vector<QuadraticIntegerBig>(
          (int)bWithMultiplicities + 1,
          QuadraticIntegerBig(iCoeffs[0], iCoeffs[1]));
    } else if (d % 8 == 1) {
      array<long int, 2> iCoeffs(
          QuadraticIntegerBig::iPellEquation(4 * iPrime));

      vector<QuadraticIntegerBig> qiFactors(
          2, QuadraticIntegerBig((iCoeffs[0] - iCoeffs[1]) / 2, iCoeffs[1]));
      qiFactors[1].conjugate();

      return qiFactors;
    } else // d % 8 == 5
      return vector<QuadraticIntegerBig>({QuadraticIntegerBig(iPrime)});
  }

  throw(string("qiFactorsRationalPrime: unknown case (d=" + std::to_string(d) +
               ", p=" + std::to_string(iPrime) + ")"));
}

double QuadraticIntegerBig::to_double() const {
  mpf_class temp(a);
  if (bIsOneMod4)
    temp += (sqrtd + 1) * b / 2;
  else
    temp += sqrtd * b;

  return temp.get_d();
}

mpz_class QuadraticIntegerBig::floor() const {
  mpf_class temp(a);
  if (bIsOneMod4)
    temp += (sqrtd + 1) * b / 2;
  else
    temp += sqrtd * b;

  return (mpz_class)::floor(temp);
}

bool QuadraticIntegerBig::isEqualTo(const int &n) const {
  return (a == n && b == 0);
}

bool QuadraticIntegerBig::isEqualTo(const AlgebraicInteger &ai) const {
  QuadraticIntegerBig qi(dynamic_cast<const QuadraticIntegerBig &>(ai));

  return (a == qi.a && b == qi.b);
}

bool QuadraticIntegerBig::isGreaterThan(const int &n) const {
  mpz_class aTemp(bIsOneMod4 ? (mpz_class(n - a) * 2 - b) : mpz_class(n - a));
  mpz_class bTemp(-b);

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

bool QuadraticIntegerBig::bIsGreaterThan(const long int &n) const {
  mpz_class aTemp(bIsOneMod4 ? (mpz_class(n - a) * 2 - b) : mpz_class(n - a));
  mpz_class bTemp(-b);

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

bool QuadraticIntegerBig::isGreaterOEThan(const int &n) const {
  return ((b == 0 && a == n) || isGreaterThan(n));
}

bool operator<(const QuadraticIntegerBig &qi1, const QuadraticIntegerBig &qi2) {
  return qi1.isLessThan(qi2);
}

bool QuadraticIntegerBig::isLessThan(const AlgebraicInteger &ai) const {
  mpz_class x(dynamic_cast<const QuadraticIntegerBig &>(ai).a);
  mpz_class y(dynamic_cast<const QuadraticIntegerBig &>(ai).b);

  mpz_class aTemp(bIsOneMod4 ? (mpz_class(a - x) * 2 + b - y)
                             : mpz_class(a - x));
  mpz_class bTemp(b - y);

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

bool QuadraticIntegerBig::isLessOEThan(const AlgebraicInteger &ai) const {
  mpz_class x(dynamic_cast<const QuadraticIntegerBig &>(ai).a);
  mpz_class y(dynamic_cast<const QuadraticIntegerBig &>(ai).b);

  if (x == a && y == b)
    return true;

  mpz_class aTemp(bIsOneMod4 ? (mpz_class(a - x) * 2 + b - y)
                             : mpz_class(a - x));
  mpz_class bTemp(b - y);

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

bool QuadraticIntegerBig::isLessThan(const int &n) const {
  mpz_class aTemp(bIsOneMod4 ? (mpz_class(a - n) * 2 + b) : mpz_class(a - n));

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

bool QuadraticIntegerBig::bIsLessThan(const long int &n) const {
  mpz_class aTemp(bIsOneMod4 ? (mpz_class(a - n) * 2 + b) : mpz_class(a - n));

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

void QuadraticIntegerBig::opp() {
  a *= -1;
  b *= -1;
}

void QuadraticIntegerBig::conjugate() {
  if (bIsOneMod4)
    a += b;

  b *= -1;
}

void QuadraticIntegerBig::substract(const AlgebraicInteger *ai) {
  a -= dynamic_cast<const QuadraticIntegerBig *>(ai)->a;
  b -= dynamic_cast<const QuadraticIntegerBig *>(ai)->b;
}

void QuadraticIntegerBig::add(const AlgebraicInteger *ai) {
  a += dynamic_cast<const QuadraticIntegerBig *>(ai)->a;
  b += dynamic_cast<const QuadraticIntegerBig *>(ai)->b;
}

void QuadraticIntegerBig::multiplyBy(const int &n) {
  a *= n;
  b *= n;
}

void QuadraticIntegerBig::multiplyBy(const AlgebraicInteger *ai) {
  mpz_class iA1(a), iB1(b);
  mpz_class iA2(dynamic_cast<const QuadraticIntegerBig *>(ai)->a);
  mpz_class iB2(dynamic_cast<const QuadraticIntegerBig *>(ai)->b);

  if (bIsOneMod4) {
    a = iA1 * iA2 + iB1 * iB2 * (d - 1) / 4;
    b = iB1 * iA2 + iA1 * iB2 + iB1 * iB2;
  } else {
    a = iA1 * iA2 + d * iB1 * iB2;
    b = iA1 * iB2 + iB1 * iA2;
  }
}

bool QuadraticIntegerBig::divideByIfDivisible(const AlgebraicInteger *ai) {
  const QuadraticIntegerBig *qi(dynamic_cast<const QuadraticIntegerBig *>(ai));

  mpz_class iA1(a), iB1(b);
  mpz_class iA2(dynamic_cast<const QuadraticIntegerBig *>(ai)->a);
  mpz_class iB2(dynamic_cast<const QuadraticIntegerBig *>(ai)->b);

  mpz_class iNorm2(qi->iNorm());

  if (this->iNorm() % iNorm2 != 0)
    return false;

  if (bIsOneMod4) {
    if ((a * (2 * iA2 + iB2) + b * iA2 - b * iB2 * (d - 1) / 2) % iNorm2 != 0)
      return false;

    a = (iA1 * (iA2 + iB2) - iB1 * iB2 * (d - 1) / 4) / iNorm2;
    b = (iB1 * iA2 - iA1 * iB2) / iNorm2;
  } else {
    if (((a * iA2 - b * iB2 * d) * 2) % iNorm2 != 0)
      return false;

    a = (iA1 * iA2 - iB1 * iB2 * d) / iNorm2;
    b = (iB1 * iA2 - iA1 * iB2) / iNorm2;
  }

  return true;
}

void QuadraticIntegerBig::divideBy(const AlgebraicInteger *ai) {
  mpz_class iA1(a), iB1(b);
  mpz_class iA2(dynamic_cast<const QuadraticIntegerBig *>(ai)->a);
  mpz_class iB2(dynamic_cast<const QuadraticIntegerBig *>(ai)->b);
  mpz_class iNorm2(dynamic_cast<const QuadraticIntegerBig *>(ai)->iNorm());

  if (bIsOneMod4) {
    a = (iA1 * (iA2 + iB2) - iB1 * iB2 * (d - 1) / 4) / iNorm2;
    b = (iB1 * iA2 - iA1 * iB2) / iNorm2;
  } else {
    a = (iA1 * iA2 - iB1 * iB2 * d) / iNorm2;
    b = (iB1 * iA2 - iA1 * iB2) / iNorm2;
  }
}

bool QuadraticIntegerBig::isDivisibleBy(const AlgebraicInteger *ai) const {
  const QuadraticIntegerBig *qi(dynamic_cast<const QuadraticIntegerBig *>(ai));

  mpz_class iNorm2(qi->iNorm());

  if (this->iNorm() % iNorm2 != 0)
    return false;

  if (bIsOneMod4) {
    if ((a * (2 * qi->a + qi->b) + b * qi->a - b * qi->b * (d - 1) / 2) %
            iNorm2 !=
        0)
      return false;
  } else {
    if (((a * qi->a - b * qi->b * d) * 2) % iNorm2 != 0)
      return false;
  }

  return true;
}

bool QuadraticIntegerBig::bIsAssociateTo(QuadraticIntegerBig qi2) {
  if (qi2.divideByIfDivisible(this)) {
    if (qi2.isInvertible())
      return true;
  }

  return false;
}

void QuadraticIntegerBig::gcd(const AlgebraicInteger *ai) {
  const QuadraticIntegerBig *qi(dynamic_cast<const QuadraticIntegerBig *>(ai));

  // -----------------------------------------
  // If one of the two is zero
  if (ai->isEqualTo(0))
    return;

  if (this->isEqualTo(0)) {
    const QuadraticIntegerBig *qi(
        dynamic_cast<const QuadraticIntegerBig *>(ai));
    a = qi->a;
    b = qi->b;
    return;
  }

  // -----------------------------------------
  // If one is bIsInvertible
  mpz_class iNorm1(abs(this->iNorm())), iNorm2(abs(qi->iNorm()));
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
  mpz_class iGCD;
  mpz_gcd(iGCD.get_mpz_t(), iNorm1.get_mpz_t(), iNorm2.get_mpz_t());
  if (iGCD == 1) {
    a = 1;
    b = 0;
    return;
  }

  vector<QuadraticIntegerBig> qiGCD(QuadraticIntegerBig(iGCD).qiPrimeFactors());
  QuadraticIntegerBig qi1(*this), qi2(*qi);

  a = 1;
  b = 0;

  for (auto qiFactor : qiGCD) {
    while (qi1.divideByIfDivisible(&qiFactor) &&
           qi2.divideByIfDivisible(&qiFactor))
      this->multiplyBy(&qiFactor);
  }
}

bool QuadraticIntegerBig::isSquareOfIvertible() const {
  if (this->isInvertible()) {
    QuadraticIntegerBig qiUnit(QuadraticIntegerBig::iFundamentalUnits[5][0],
                               QuadraticIntegerBig::iFundamentalUnits[5][1]);
    QuadraticIntegerBig qiTemp(*this);
    unsigned int iPower(0);

    while (!qiTemp.isEqualTo(1)) {
      qiTemp.divideBy(&qiUnit);
      iPower++;
    }

    return (iPower % 2 == 0);
  }
  return false;
}

bool QuadraticIntegerBig::isInvertible() const {
  mpz_class iNorm(this->iNorm());
  return (iNorm == 1 || iNorm == -1);
}

void QuadraticIntegerBig::removeSquareFactors() {
  map<QuadraticIntegerBig, unsigned int> qiDecomp(this->qiPrimeDecomposition());
  for (auto it : qiDecomp) {
    unsigned int iMax(it.second % 2 ? it.second - 1 : it.second);
    for (unsigned int i(1); i <= iMax; i++)
      this->divideBy(&it.first);
  }
}

void QuadraticIntegerBig::set(AlgebraicInteger *ai) {
  QuadraticIntegerBig *qi(dynamic_cast<QuadraticIntegerBig *>(ai));

  a = qi->a;
  b = qi->b;
}

void QuadraticIntegerBig::set(const int &n) {
  a = n;
  b = 0;
}

mpz_class QuadraticIntegerBig::iNorm() const {
  if (bIsOneMod4)
    return (a * (a + b) - b * b * (d - 1) / 4);
  else
    return (a * a - b * b * mpz_class(d));
}

mpz_class QuadraticIntegerBig::iTrace() const {
  if (bIsOneMod4)
    return (a * 2 + b);
  else
    return a * 2;
}

mpz_class
QuadraticIntegerBig::iSQRTsup_quotient(const QuadraticIntegerBig &qiNum,
                                       const QuadraticIntegerBig &qiDen) {
  mpz_class x, y, z, s;
  mpz_class iNumTrace(qiNum.iTrace()), iDenTrace(qiDen.iTrace());

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
    mpz_class is;
    if (y > 0) {
      is = y * y * d - 1;
      mpz_sqrt(is.get_mpz_t(), is.get_mpz_t());
      is++; // Now, is = iSQRTsup( y * y * d )
    } else
      is = 0;

    s = sqrtSupQuotient(x + is, z);

    while ((z * (s + 1) * (s + 1) - x) * (z * (s + 1) * (s + 1) - x) <=
           y * y * d)
      s--;
  } else {
    mpz_class is(y * y * d);
    mpz_sqrt(is.get_mpz_t(), is.get_mpz_t());

    s = sqrtSupQuotient(x - is, z);

    mpz_class iTemp(x - (s - 1) * (s - 1) * z);

    while (iTemp < 0 || iTemp * iTemp <= is) {
      s--;
      iTemp = x - (s - 1) * (s - 1) * z;
    }
  }

  return s;
}

mpz_class
QuadraticIntegerBig::iSQRT_quotient(const QuadraticIntegerBig &qiNum,
                                    const QuadraticIntegerBig &qiDen) {
  mpz_class x, y, z, s;
  mpz_class iNumTrace(qiNum.iTrace()), iDenTrace(qiDen.iTrace());

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

  mpz_class is(y * y * d);
  mpz_sqrt(is.get_mpz_t(), is.get_mpz_t());
  if (y >= 0) {
    s = sqrtQuotient(x + is, z);

    while ((z * (s + 1) * (s + 1) - x) * (z * (s + 1) * (s + 1) - x) <=
           y * y * d)
      s++;
  } else {
    s = sqrtSupQuotient(x - is, z);

    mpz_class iTemp(z * s * s - x);

    while ((iTemp > 0) || (iTemp <= 0 && iTemp * iTemp < y * y * d)) {
      s--;
      iTemp = z * s * s - x;
    }
  }

  return s;
}

AlgebraicInteger *QuadraticIntegerBig::copy() const {
  return new QuadraticIntegerBig(*this);
}

AlgebraicInteger *QuadraticIntegerBig::aiCopyToInteger(const int &n) const {
  return new QuadraticIntegerBig(n);
}

std::ostream &QuadraticIntegerBig::print(std::ostream &o) const {
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

std::string QuadraticIntegerBig::to_string(const string &strFormat,
                                           const bool &bProtect) const {
  mpz_class babs(abs(b));

  if (strFormat == "mathematica") {
    if (bIsOneMod4) {
      if (a == 0) {
        if (b == 0)
          return "0";
        else
          return b.get_str() + "*(1+Sqrt[" + std::to_string(d) + "])/2";
      } else {
        if (b == 0)
          return a.get_str();
        else
          return (bProtect ? "(" : "") + a.get_str() + (b < 0 ? "-" : "+") +
                 babs.get_str() + "*(1+Sqrt[" + std::to_string(d) +
                 (bProtect ? "])/2)" : "])/2");
      }
    } else {
      if (a == 0) {
        if (b == 0)
          return "0";
        else
          return b.get_str() + "*Sqrt[" + std::to_string(d) + "]";
      } else {
        if (b == 0)
          return a.get_str();
        else
          return (bProtect ? "(" : "") + a.get_str() + (b < 0 ? "-" : "+") +
                 babs.get_str() + "*Sqrt[" + std::to_string(d) +
                 (bProtect ? "])" : "]");
      }
    }

  } else if (strFormat == "filename") {
    if (a == 0) {
      if (b == 0)
        return "0";
      else
        return (b == 1 ? "" : b.get_str()) + "T" + std::to_string(d);
    } else {
      if (b == 0)
        return a.get_str();
      else
        return (bProtect ? "(" : "") + a.get_str() + (b < 0 ? "-" : "+") +
               (abs(b) == 1 ? "" : babs.get_str()) + "T" + std::to_string(d) +
               (bProtect ? ")" : "");
    }
  } else {
    if (a == 0) {
      if (b == 0)
        return "0";
      else
        return b.get_str() + "*T" + std::to_string(d);
    } else {
      if (b == 0)
        return a.get_str();
      else
        return (bProtect ? "(" : "") + a.get_str() + (b < 0 ? "-" : "+") +
               babs.get_str() + "*T" + std::to_string(d) +
               (bProtect ? ")" : "");
    }
  }
}

string QuadraticIntegerBig::get_classname() const {
  return "QuadraticIntegerBig";
}

// ------------------------------------------------------------------------
// Operators
QuadraticIntegerBig &
QuadraticIntegerBig::operator=(const QuadraticIntegerBig &qi) {
  a = qi.a;
  b = qi.b;

  return *this;
}

QuadraticIntegerBig &
QuadraticIntegerBig::operator/=(const QuadraticIntegerBig &qi) {
  this->divideBy(&qi);
  return *this;
}

QuadraticIntegerBig &
QuadraticIntegerBig::operator*=(const QuadraticIntegerBig &qi) {
  this->multiplyBy(&qi);
  return *this;
}

QuadraticIntegerBig
QuadraticIntegerBig::operator+(const QuadraticIntegerBig &qi) const {
  return QuadraticIntegerBig(a + qi.a, b + qi.b);
}

QuadraticIntegerBig
QuadraticIntegerBig::operator*(const QuadraticIntegerBig &qi) const {
  QuadraticIntegerBig qiRes(*this);
  qiRes.multiplyBy(&qi);

  return qiRes;
}

QuadraticIntegerBig
QuadraticIntegerBig::operator-(const QuadraticIntegerBig &qi) const {
  return QuadraticIntegerBig(a - qi.a, b - qi.b);
}

QuadraticIntegerBig QuadraticIntegerBig::operator-() const {
  return QuadraticIntegerBig(-a, -b);
}

bool QuadraticIntegerBig::operator>(const QuadraticIntegerBig &ri) const {
  return !isLessOEThan(ri);
}

bool QuadraticIntegerBig::operator==(const QuadraticIntegerBig &qi) const {
  return (a == qi.a && b == qi.b);
}

bool QuadraticIntegerBig::operator==(const long int &n) const {
  return (a == n && b == 0);
}
