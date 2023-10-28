#include "rcyclotomic7integer.h"

using namespace RCyclotomic7Integer_constants;

vector<vector<int>> RCyclotomic7Integer::iPowersL = vector<vector<int>>(0);

array<array<array<mpq_class, 2>, 3>, 2>
    RCyclotomic7Integer::mpqRationalApproximations = {
        {{{{{mpq_class(mpz_class("21689544433045364711933845"),
                       mpz_class("17393664153274833967747621")),
             mpq_class(mpz_class("10011913795355660647342091"),
                       mpz_class("8028931480120743221914558"))}},
           {{mpq_class(mpz_class("-17393664153274833967747621"),
                       mpz_class("39083208586320198679681466")),
             mpq_class(mpz_class("-6693130287087395697996053"),
                       mpz_class("15039327240109012928088481"))}},
           {{mpq_class(mpz_class("-39083208586320198679681466"),
                       mpz_class("21689544433045364711933845")),
             mpq_class(mpz_class("-18040845275476403869256649"),
                       mpz_class("10011913795355660647342091"))}}}},

         {{{{mpq_class(
                 mpz_class(
                     "118617165482138066506657879508313239222975636291488524417"
                     "50139798979225553949651060857197331644225179"),
                 mpz_class(
                     "951235811143336161914788967901403952137597001530575172344"
                     "6393067643895366950437131247228319682746513")),
             mpq_class(
                 mpz_class(
                     "4279766040305502977765390182184089005892462452168265202757047478106033\
2669976203083482085854272332411"),
                 mpz_class(
                     "343210588813703318689430030482808282837641326940761405007"
                     "09644855433530765149420395818437777189360798"))}},
           {{mpq_class(
                 mpz_class(
                     "-95123581114333616191478896790140395213759700153057517234"
                     "46393067643895366950437131247228319682746513"),
                 mpz_class(
                     "2137407465964716826981367762984536344367353364445460416519653286662312\
0920900088192104425651326971692")),
             mpq_class(
                 mpz_class(
                     "-36001241779569731130349525376494178697411246448174740477"
                     "11165575344850787496486426221461372761653920761"),
                 mpz_class(
                     "8089405598719431309982823636977886334493376006819893508897776026292217\
217866533715958193087866151300015"))}},
           {{mpq_class(
                 mpz_class(
                     "-21374074659647168269813677629845363443673533644454604165"
                     "196532866623120920900088192104425651326971692"),
                 mpz_class(
                     "1186171654821380665066578795083132392229756362914885244175013979897922\
5553949651060857197331644225179")),
             mpq_class(
                 mpz_class(
                     "-77118719284425361646596904870121718342688757215758792528"
                     "280119636493863435125623479300523631461693209"),
                 mpz_class(
                     "427976604030550297776539018218408900589246245216826520275"
                     "70474781060332669976203083482085854272332411"))}}}}}};

bool RCyclotomic7Integer::bClassInitialized = false;

void RCyclotomic7Integer::initialize() {
  if (bClassInitialized)
    return;

  iPowersL = vector<vector<int>>({{-2, -1, -2}, {3, 0, 1}});

  // pari_init(1000000,2); TODO

  bClassInitialized = true;
}

void RCyclotomic7Integer::computePrimesDecomposition_factorMinimalPolynomial(
    GEN &gFactors, const unsigned int &iPrime) {
  vector<int> iMinPol({-1, -2, 1, 1});
  unsigned int iMinimalPolynomialSize(4);

  GEN gPrime(stoi(iPrime));

  const long iVar = 0; // Polynomial variable
  long iCoefficient;

  // ----------------------------------------------
  // Create the polynomial mod iPrime
  GEN pol(cgetg(
      iMinimalPolynomialSize + 2,
      t_POL)); // iMinimalPolynomialSize coefficients, sign and variable number
  pol[1] = evalsigne(1) | evalvarn(iVar) |
           evallgefint(iMinimalPolynomialSize +
                       2); // not equal to zero, variable, coefficients

  for (unsigned int i(0); i < iMinimalPolynomialSize; i++) {
    iCoefficient = iMinPol[i];

    if (iCoefficient <
        0) // We add some multiple of iPrime so that iCoefficient is positive
      iCoefficient +=
          iPrime * ((-iCoefficient % iPrime) ? -iCoefficient / iPrime + 1
                                             : -iCoefficient / iPrime);

    pol[i + 2] = (long)stoi(iCoefficient % iPrime);
  }

  pol = RgX_to_FpX(pol, gPrime);      // Reduction mod gPrime
  gFactors = FpX_factor(pol, gPrime); // Facorization
}

bool RCyclotomic7Integer::computePrimeDecomposition(
    const unsigned int &iPrime) {
  pari_sp av = avma; // current state of the PARI stack

  // -----------------------------------------------
  // Factorization of the minimal polynomial mod F_p
  GEN gFactors, gExponents;
  computePrimesDecomposition_factorMinimalPolynomial(gFactors, iPrime);

  gExponents = gel(gFactors, 2); // Exponents of the irreducible factors
  gFactors = gel(gFactors, 1);   // Irreducible factors
  unsigned int iFactorsCount(
      (unsigned int)(lg(gFactors) - 1)); // Number of irreducible factors
  // The above are not really useful since gExponents * gFactors * iFactorsCount
  // = 3.

  if (iFactorsCount == 1 && gExponents[1] == 1) // iPrime is prime in O_K
  {
    iPrimesDecomposition[iPrime] =
        vector<RCyclotomic7Integer>({RCyclotomic7Integer(iPrime)});
  } else {
    iPrimesDecomposition[iPrime] = vector<RCyclotomic7Integer>(0);
    long int iC0, iC1, iNorm, iResidualNorm;
    for (unsigned int i(1); i <= iFactorsCount &&
                            iPrimesDecomposition[iPrime].size() < iFactorsCount;
         i++) {
      // The polynomial is of degree 1
      iC1 = gtolong(gel(gel(gFactors, i), 3));

      for (int k(-2); k <= 1; k++) {
        iC0 = gtolong(gel(gel(gFactors, i), 2)) + k * (long int)iPrime;
        RCyclotomic7Integer rci(-iC0 + iC1, -iC0, -iC0);

        mpz_class mpz_norm(rci.iNorm());
        if (!mpz_norm.fits_slong_p())
          throw(string("RCyclotomic7Integer::computePrimeDecomposition: Number "
                       "too big"));
        iNorm = abs(mpz_norm.get_si());

        iResidualNorm = iNorm / ugcd(iPrime, iNorm);

        if (iResidualNorm == 1 || iResidualNorm == -1)
          computePrimesDecomposition_addFactor(iPrime, iFactorsCount, rci);

        if (iFactorsCount == iPrimesDecomposition[iPrime].size())
          break;

        vector<unsigned long> iResidualPrimeFactors(
            primeFactors<unsigned long>(iResidualNorm));

        bool bResidualFactorsFactorized(true);
        for (auto iRes : iResidualPrimeFactors) {
          if (iPrimesDecomposition.find(iRes) == iPrimesDecomposition.end()) {
            bResidualFactorsFactorized = false;
            break;
          }
        }

        if (!bResidualFactorsFactorized)
          continue;

        for (auto iPFactor : iResidualPrimeFactors) {
          for (auto iF : iPrimesDecomposition[iPFactor]) {
            while (rci.divideByIfDivisible(&iF))
              ;
          }
        }

        iNorm = abs(rci.iNorm().get_si());
        iResidualNorm = iNorm / ugcd(iPrime, iNorm);

        if (iResidualNorm == 1)
          computePrimesDecomposition_addFactor(iPrime, iFactorsCount, rci);

        if (iFactorsCount == iPrimesDecomposition[iPrime].size())
          break;
      }
    }
  }

  avma = av; // cleaning the stack

  if (iFactorsCount != iPrimesDecomposition[iPrime].size()) {
    iPrimesDecomposition[iPrime].clear();
    return false;
  }

  return true;
}

void RCyclotomic7Integer::computePrimesDecomposition_addFactor(
    const unsigned int &iPrime, const unsigned int &iFactorsCount,
    RCyclotomic7Integer rci) {
  vector<RCyclotomic7Integer> *ptrList(&iPrimesDecomposition[iPrime]);
  unsigned int iListSize(ptrList->size()), i;

  for (unsigned int j(0); j < 3; j++) {
    // We already have all the factors
    if (iListSize >= iFactorsCount)
      return;

    for (i = 0; i < iListSize; i++) {
      if ((*ptrList)[i].isAssociateTo(&rci))
        break;
    }

    if (i == iListSize) {
      if (rci.isLessThan(0))
        rci.multiplyBy(-1);

      ptrList->push_back(rci);
      iListSize++;
    }

    rci.conjugate(2);
  }
}

RCyclotomic7Integer::RCyclotomic7Integer() {
  if (!bClassInitialized)
    initialize();

  iC = array<mpz_class, 3>({0, 0, 0});
}

RCyclotomic7Integer::RCyclotomic7Integer(const RCyclotomic7Integer &rci) {
  iC = rci.iC;
}

RCyclotomic7Integer::RCyclotomic7Integer(const int &n) {
  iC = array<mpz_class, 3>({-n, -n, -n});
}

RCyclotomic7Integer::RCyclotomic7Integer(const int &c1, const int &c2, const int &c3) {
  iC = array<mpz_class, 3>({c1, c2, c3});
}


RCyclotomic7Integer::RCyclotomic7Integer(const mpz_class &n) {
  iC = array<mpz_class, 3>({-n, -n, -n});
}

RCyclotomic7Integer::RCyclotomic7Integer(
    const array<mpz_class, 3> &iCoefficients)
    : iC(iCoefficients) {}

RCyclotomic7Integer::~RCyclotomic7Integer() {}

std::string RCyclotomic7Integer::to_string(const std::string &strFormat,
                                           const bool &bProtect) const {
  if (iC[0] == iC[1] && iC[1] == iC[2]) {
    mpz_class i(-iC[0]);
    return i.get_str();
  } else {
    if (strFormat == "mathematica")
      return string("RCI7[{" + iC[0].get_str() + "," + iC[1].get_str() + "," +
                    iC[2].get_str() + "}]");
    else if (strFormat == "pari")
      return string("RCI7(" + iC[0].get_str() + "," + iC[1].get_str() + "," +
                    iC[2].get_str() + ")");
    else {
      string str("[");
      for (unsigned int i(0); i < 3; i++)
        str += (i ? "," : "") + iC[i].get_str();
      str += "]";
      return str;
    }
  }
}

double RCyclotomic7Integer::to_double() const {
  return (dl1 * iC[0].get_d() + dl2 * iC[1].get_d() + dl3 * iC[2].get_d());
}

interval RCyclotomic7Integer::to_interval() const {
  if (iC[0] == iC[1] && iC[1] == iC[2])
    return interval(gaol_intervalFromMPZclass(-iC[0]));

  return (interval(gaol_intervalFromMPZclass(iC[0])) * gaol_l1 +
          interval(gaol_intervalFromMPZclass(iC[1])) * gaol_l2 +
          interval(gaol_intervalFromMPZclass(iC[2])) * gaol_l3);
}

long int RCyclotomic7Integer::floor() const {
  interval gaol(gaol::floor(to_interval()));

  if (!gaol.is_an_int()) {
    // TODO
    cout << "Gaol: " << gaol << "; " << *this << endl;
    throw(string("RCyclotomic7Integer::to_interval() : Not an int"));
  }

  return gaol.left();
}

mpz_class RCyclotomic7Integer::iNorm() const {
  return (iC[0] * iC[0] * iC[0] + iC[1] * iC[1] * iC[1] +
          iC[2] * iC[2] * iC[2] - iC[0] * iC[1] * iC[2] +
          3 * (iC[0] * iC[0] * iC[1] + iC[1] * iC[1] * iC[2] +
               iC[0] * iC[2] * iC[2]) -
          4 * (iC[0] * iC[1] * iC[1] + iC[0] * iC[0] * iC[2] +
               iC[1] * iC[2] * iC[2]));
}

bool RCyclotomic7Integer::isLongInt() const {
  if (iC[0] == iC[1] && iC[1] == iC[2])
    return iC[0].fits_slong_p();

  return false;
}

void RCyclotomic7Integer::conjugate(unsigned int i) {
  i %= 3;

  mpz_class a0(iC[0]), a1(iC[1]), a2(iC[2]);

  if (i == 2) // second element of the Galois group
  {
    iC[0] = a2;
    iC[1] = a0;
    iC[2] = a1;
  } else if (i == 0) // third element of the Galois group
  {
    iC[0] = a1;
    iC[1] = a2;
    iC[2] = a0;
  }
  // else: id
}

std::ostream &RCyclotomic7Integer::print(std::ostream &o) const {
  if (iC[0] == iC[1] && iC[1] == iC[2])
    o << (-1 * iC[0]);
  else {
    o << "[";
    for (unsigned int i(0); i < 3; i++)
      o << (i ? "," : "") << iC[i];
    o << "]";
  }

  return o;
}

bool RCyclotomic7Integer::isEqualTo(const int &n) const {
  return (iC[0] == -n && iC[1] == -n && iC[2] == -n);
}

bool RCyclotomic7Integer::isEqualTo(const AlgebraicInteger &ai) const {
  return iC == dynamic_cast<const RCyclotomic7Integer &>(ai).iC;
}

bool RCyclotomic7Integer::isGreaterOEThan(const int &n) const {
  if (iC[0] == -n && iC[1] == -n && iC[2] == -n)
    return true;

  interval i(interval(gaol_intervalFromMPZclass(iC[0])) * gaol_l1 +
             interval(gaol_intervalFromMPZclass(iC[1])) * gaol_l2 +
             interval(gaol_intervalFromMPZclass(iC[2])) * gaol_l3 -
             interval(n));

  if (i.certainly_negative())
    return false;

  if (i.certainly_strictly_positive())
    return true;

  throw(string("RCyclotomic7Integer::isGreaterOEThan : Cannot decide"));
}

bool RCyclotomic7Integer::isGreaterThan(const int &n) const {
  interval i(interval(gaol_intervalFromMPZclass(iC[0])) * gaol_l1 +
             interval(gaol_intervalFromMPZclass(iC[1])) * gaol_l2 +
             interval(gaol_intervalFromMPZclass(iC[2])) * gaol_l3 -
             interval(n));

  if (i.certainly_negative())
    return false;

  if (i.certainly_strictly_positive())
    return true;

  throw(string("RCyclotomic7Integer::isGreaterThan : Cannot decide"));
}

bool RCyclotomic7Integer::isLessOEThan(const AlgebraicInteger &ai) const {
  RCyclotomic7Integer rci(dynamic_cast<const RCyclotomic7Integer &>(ai));

  if (iC == rci.iC)
    return true;

  interval i(interval(gaol_intervalFromMPZclass(iC[0] - rci.iC[0])) * gaol_l1 +
             interval(gaol_intervalFromMPZclass(iC[1] - rci.iC[1])) * gaol_l2 +
             interval(gaol_intervalFromMPZclass(iC[2] - rci.iC[2])) * gaol_l3);

  if (i.certainly_positive())
    return false;

  if (i.certainly_strictly_negative())
    return true;

  mpz_class i0(iC[0] - rci.iC[0]), i1(iC[1] - rci.iC[1]), i2(iC[2] - rci.iC[2]);

  // ---------------------------------------------------
  // Here, we will try rational approximation (10^-50)
  mpq_class mpqInf((i0 > 0 ? mpqRationalApproximations[0][0][0]
                           : mpqRationalApproximations[0][0][1]) *
                   i0);
  mpqInf += (i1 > 0 ? mpqRationalApproximations[0][1][0]
                    : mpqRationalApproximations[0][1][1]) *
            i1;
  mpqInf += (i2 > 0 ? mpqRationalApproximations[0][2][0]
                    : mpqRationalApproximations[0][2][1]) *
            i2;

  if (mpqInf > 0)
    return false;

  mpq_class mpqSup((i0 > 0 ? mpqRationalApproximations[0][0][1]
                           : mpqRationalApproximations[0][0][0]) *
                   i0);
  mpqSup += (i1 > 0 ? mpqRationalApproximations[0][1][1]
                    : mpqRationalApproximations[0][1][0]) *
            i1;
  mpqSup += (i2 > 0 ? mpqRationalApproximations[0][2][1]
                    : mpqRationalApproximations[0][2][0]) *
            i2;

  if (mpqSup < 0)
    return true;

  // ---------------------------------------------------
  // Here, we will try rational approximation (10^-200)
  mpqInf = (i0 > 0 ? mpqRationalApproximations[1][0][0]
                   : mpqRationalApproximations[1][0][1]) *
           i0;
  mpqInf += (i1 > 0 ? mpqRationalApproximations[1][1][0]
                    : mpqRationalApproximations[1][1][1]) *
            i1;
  mpqInf += (i2 > 0 ? mpqRationalApproximations[1][2][0]
                    : mpqRationalApproximations[1][2][1]) *
            i2;

  if (mpqInf > 0)
    return false;

  mpqSup = (i0 > 0 ? mpqRationalApproximations[1][0][1]
                   : mpqRationalApproximations[1][0][0]) *
           i0;
  mpqSup += (i1 > 0 ? mpqRationalApproximations[1][1][1]
                    : mpqRationalApproximations[1][1][0]) *
            i1;
  mpqSup += (i2 > 0 ? mpqRationalApproximations[1][2][1]
                    : mpqRationalApproximations[1][2][0]) *
            i2;

  if (mpqSup < 0)
    return true;

  throw(string("RCyclotomic7Integer::isLessOEThan : Cannot decide"));
}

bool operator<(const RCyclotomic7Integer &rci1,
               const RCyclotomic7Integer &rci2) {
  return rci1.isLessThan(rci2);
}

bool RCyclotomic7Integer::isLessThan(const AlgebraicInteger &ai) const {
  RCyclotomic7Integer rci(dynamic_cast<const RCyclotomic7Integer &>(ai));

  interval i(interval(gaol_intervalFromMPZclass(iC[0] - rci.iC[0])) * gaol_l1 +
             interval(gaol_intervalFromMPZclass(iC[1] - rci.iC[1])) * gaol_l2 +
             interval(gaol_intervalFromMPZclass(iC[2] - rci.iC[2])) * gaol_l3);

  if (i.certainly_positive())
    return false;

  if (i.certainly_strictly_negative())
    return true;

  throw(string("RCyclotomic7Integer::isLessThan(ai) : Cannot decide"));
}

bool RCyclotomic7Integer::isLessThan(const int &n) const {
  interval i(interval(gaol_intervalFromMPZclass(iC[0])) * gaol_l1 +
             interval(gaol_intervalFromMPZclass(iC[1])) * gaol_l2 +
             interval(gaol_intervalFromMPZclass(iC[2])) * gaol_l3 -
             interval(n));

  if (i.certainly_positive())
    return false;

  if (i.certainly_strictly_negative())
    return true;

  mpz_class i0(iC[0] + n), i1(iC[1] + n), i2(iC[2] + n);

  // ---------------------------------------------------
  // Here, we will try rational approximation (10^-50)
  mpq_class mpqInf((i0 > 0 ? mpqRationalApproximations[0][0][0]
                           : mpqRationalApproximations[0][0][1]) *
                   i0);
  mpqInf += (i1 > 0 ? mpqRationalApproximations[0][1][0]
                    : mpqRationalApproximations[0][1][1]) *
            i1;
  mpqInf += (i2 > 0 ? mpqRationalApproximations[0][2][0]
                    : mpqRationalApproximations[0][2][1]) *
            i2;

  if (mpqInf > 0)
    return false;

  mpq_class mpqSup((i0 > 0 ? mpqRationalApproximations[0][0][1]
                           : mpqRationalApproximations[0][0][0]) *
                   i0);
  mpqSup += (i1 > 0 ? mpqRationalApproximations[0][1][1]
                    : mpqRationalApproximations[0][1][0]) *
            i1;
  mpqSup += (i2 > 0 ? mpqRationalApproximations[0][2][1]
                    : mpqRationalApproximations[0][2][0]) *
            i2;

  if (mpqSup < 0)
    return true;

  // ---------------------------------------------------
  // Here, we will try rational approximation (10^-200)
  mpqInf = (i0 > 0 ? mpqRationalApproximations[1][0][0]
                   : mpqRationalApproximations[1][0][1]) *
           i0;
  mpqInf += (i1 > 0 ? mpqRationalApproximations[1][1][0]
                    : mpqRationalApproximations[1][1][1]) *
            i1;
  mpqInf += (i2 > 0 ? mpqRationalApproximations[1][2][0]
                    : mpqRationalApproximations[1][2][1]) *
            i2;

  if (mpqInf > 0)
    return false;

  mpqSup = (i0 > 0 ? mpqRationalApproximations[1][0][1]
                   : mpqRationalApproximations[1][0][0]) *
           i0;
  mpqSup += (i1 > 0 ? mpqRationalApproximations[1][1][1]
                    : mpqRationalApproximations[1][1][0]) *
            i1;
  mpqSup += (i2 > 0 ? mpqRationalApproximations[1][2][1]
                    : mpqRationalApproximations[1][2][0]) *
            i2;

  if (mpqSup < 0)
    return true;

  throw(string("RCyclotomic7Integer::isLessThan(int) : Cannot decide"));
}

void RCyclotomic7Integer::opp() {
  iC[0] *= -1;
  iC[1] *= -1;
  iC[2] *= -1;
}

void RCyclotomic7Integer::substract(const AlgebraicInteger *ai) {
  const RCyclotomic7Integer *rci(dynamic_cast<const RCyclotomic7Integer *>(ai));

  iC[0] -= rci->iC[0];
  iC[1] -= rci->iC[1];
  iC[2] -= rci->iC[2];
}

void RCyclotomic7Integer::add(const AlgebraicInteger *ai) {
  const RCyclotomic7Integer *rci(dynamic_cast<const RCyclotomic7Integer *>(ai));

  iC[0] += rci->iC[0];
  iC[1] += rci->iC[1];
  iC[2] += rci->iC[2];
}

void RCyclotomic7Integer::multiplyBy(const AlgebraicInteger *ai) {
  const RCyclotomic7Integer *rci(dynamic_cast<const RCyclotomic7Integer *>(ai));

  mpz_class a0(iC[0]), a1(iC[1]), a2(iC[2]);
  mpz_class b0(rci->iC[0]), b1(rci->iC[1]), b2(rci->iC[2]);

  iC[0] = -2 * (a0 * b0 + a1 * b1) + a1 * b0 + (a0 + a2) * b1 + (a1 - a2) * b2;
  iC[1] = (a2 - a0) * b0 - 2 * (a1 * b1 + a2 * b2) + a2 * b1 + (a0 + a1) * b2;
  iC[2] = -2 * (a0 * b0 + a2 * b2) + (a1 + a2) * b0 + (a0 - a1) * b1 + a0 * b2;
}

void RCyclotomic7Integer::multiplyBy(const int &n) {
  iC[0] *= n;
  iC[1] *= n;
  iC[2] *= n;
}

bool RCyclotomic7Integer::divideByIfDivisible(const AlgebraicInteger *ai) {
  const RCyclotomic7Integer *rci(dynamic_cast<const RCyclotomic7Integer *>(ai));

  mpz_class a0(iC[0]), a1(iC[1]), a2(iC[2]);
  mpz_class b0(rci->iC[0]), b1(rci->iC[1]), b2(rci->iC[2]);

  mpz_class iNormb(rci->iNorm());
  // TODO: ré-écrire pour otpimiser
  if (((-a0 * b0 * b0 - a1 * b0 * b0 + a2 * b0 * b0 - 2 * a0 * b0 * b1 +
        3 * a1 * b0 * b1 - a2 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 +
        3 * a0 * b0 * b2 - a2 * b0 * b2 + 2 * a0 * b1 * b2 - 3 * a1 * b1 * b2 +
        2 * a2 * b1 * b2 - 2 * a0 * b2 * b2 + 2 * a1 * b2 * b2 - a2 * b2 * b2) %
       iNormb) != 0 ||
      ((-a0 * b0 * b0 - 2 * a1 * b0 * b0 + 2 * a2 * b0 * b0 - a0 * b0 * b1 +
        3 * a1 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 - a2 * b1 * b1 +
        2 * a0 * b0 * b2 + 2 * a1 * b0 * b2 - 3 * a2 * b0 * b2 - a0 * b1 * b2 -
        2 * a1 * b1 * b2 + 3 * a2 * b1 * b2 + a1 * b2 * b2 - a2 * b2 * b2) %
       iNormb) != 0 ||
      ((-a0 * b0 * b0 + a2 * b0 * b0 - 3 * a0 * b0 * b1 + 2 * a1 * b0 * b1 +
        2 * a2 * b0 * b1 + 2 * a0 * b1 * b1 - a1 * b1 * b1 - 2 * a2 * b1 * b1 +
        3 * a0 * b0 * b2 - a1 * b0 * b2 - 2 * a2 * b0 * b2 - a1 * b1 * b2 +
        3 * a2 * b1 * b2 - a0 * b2 * b2 + a1 * b2 * b2 - a2 * b2 * b2) %
       iNormb) != 0)
    return false;

  iC[0] =
      (-a0 * b0 * b0 - a1 * b0 * b0 + a2 * b0 * b0 - 2 * a0 * b0 * b1 +
       3 * a1 * b0 * b1 - a2 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 +
       3 * a0 * b0 * b2 - a2 * b0 * b2 + 2 * a0 * b1 * b2 - 3 * a1 * b1 * b2 +
       2 * a2 * b1 * b2 - 2 * a0 * b2 * b2 + 2 * a1 * b2 * b2 - a2 * b2 * b2) /
      iNormb;

  iC[1] =
      (-a0 * b0 * b0 - 2 * a1 * b0 * b0 + 2 * a2 * b0 * b0 - a0 * b0 * b1 +
       3 * a1 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 - a2 * b1 * b1 +
       2 * a0 * b0 * b2 + 2 * a1 * b0 * b2 - 3 * a2 * b0 * b2 - a0 * b1 * b2 -
       2 * a1 * b1 * b2 + 3 * a2 * b1 * b2 + a1 * b2 * b2 - a2 * b2 * b2) /
      iNormb;

  iC[2] =
      (-a0 * b0 * b0 + a2 * b0 * b0 - 3 * a0 * b0 * b1 + 2 * a1 * b0 * b1 +
       2 * a2 * b0 * b1 + 2 * a0 * b1 * b1 - a1 * b1 * b1 - 2 * a2 * b1 * b1 +
       3 * a0 * b0 * b2 - a1 * b0 * b2 - 2 * a2 * b0 * b2 - a1 * b1 * b2 +
       3 * a2 * b1 * b2 - a0 * b2 * b2 + a1 * b2 * b2 - a2 * b2 * b2) /
      iNormb;

  return true;
}

void RCyclotomic7Integer::divideBy(const AlgebraicInteger *ai) {
  const RCyclotomic7Integer *rci(dynamic_cast<const RCyclotomic7Integer *>(ai));

  mpz_class a0(iC[0]), a1(iC[1]), a2(iC[2]);
  mpz_class b0(rci->iC[0]), b1(rci->iC[1]), b2(rci->iC[2]);

  mpz_class iNormb(rci->iNorm());

  iC[0] =
      (-a0 * b0 * b0 - a1 * b0 * b0 + a2 * b0 * b0 - 2 * a0 * b0 * b1 +
       3 * a1 * b0 * b1 - a2 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 +
       3 * a0 * b0 * b2 - a2 * b0 * b2 + 2 * a0 * b1 * b2 - 3 * a1 * b1 * b2 +
       2 * a2 * b1 * b2 - 2 * a0 * b2 * b2 + 2 * a1 * b2 * b2 - a2 * b2 * b2) /
      iNormb;

  iC[1] =
      (-a0 * b0 * b0 - 2 * a1 * b0 * b0 + 2 * a2 * b0 * b0 - a0 * b0 * b1 +
       3 * a1 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 - a2 * b1 * b1 +
       2 * a0 * b0 * b2 + 2 * a1 * b0 * b2 - 3 * a2 * b0 * b2 - a0 * b1 * b2 -
       2 * a1 * b1 * b2 + 3 * a2 * b1 * b2 + a1 * b2 * b2 - a2 * b2 * b2) /
      iNormb;

  iC[2] =
      (-a0 * b0 * b0 + a2 * b0 * b0 - 3 * a0 * b0 * b1 + 2 * a1 * b0 * b1 +
       2 * a2 * b0 * b1 + 2 * a0 * b1 * b1 - a1 * b1 * b1 - 2 * a2 * b1 * b1 +
       3 * a0 * b0 * b2 - a1 * b0 * b2 - 2 * a2 * b0 * b2 - a1 * b1 * b2 +
       3 * a2 * b1 * b2 - a0 * b2 * b2 + a1 * b2 * b2 - a2 * b2 * b2) /
      iNormb;
}

bool RCyclotomic7Integer::isAssociateTo(const AlgebraicInteger *ai) {
  const RCyclotomic7Integer *rci(dynamic_cast<const RCyclotomic7Integer *>(ai));
  mpz_class iNormb(rci->iNorm());

  if (iNorm() != iNormb)
    return false;

  mpz_class a0(iC[0]), a1(iC[1]), a2(iC[2]);
  mpz_class b0(rci->iC[0]), b1(rci->iC[1]), b2(rci->iC[2]);

  return (
      ((-a0 * b0 * b0 - a1 * b0 * b0 + a2 * b0 * b0 - 2 * a0 * b0 * b1 +
        3 * a1 * b0 * b1 - a2 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 +
        3 * a0 * b0 * b2 - a2 * b0 * b2 + 2 * a0 * b1 * b2 - 3 * a1 * b1 * b2 +
        2 * a2 * b1 * b2 - 2 * a0 * b2 * b2 + 2 * a1 * b2 * b2 - a2 * b2 * b2) %
       iNormb) == 0 &&
      ((-a0 * b0 * b0 - 2 * a1 * b0 * b0 + 2 * a2 * b0 * b0 - a0 * b0 * b1 +
        3 * a1 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 - a2 * b1 * b1 +
        2 * a0 * b0 * b2 + 2 * a1 * b0 * b2 - 3 * a2 * b0 * b2 - a0 * b1 * b2 -
        2 * a1 * b1 * b2 + 3 * a2 * b1 * b2 + a1 * b2 * b2 - a2 * b2 * b2) %
       iNormb) == 0 &&
      ((-a0 * b0 * b0 + a2 * b0 * b0 - 3 * a0 * b0 * b1 + 2 * a1 * b0 * b1 +
        2 * a2 * b0 * b1 + 2 * a0 * b1 * b1 - a1 * b1 * b1 - 2 * a2 * b1 * b1 +
        3 * a0 * b0 * b2 - a1 * b0 * b2 - 2 * a2 * b0 * b2 - a1 * b1 * b2 +
        3 * a2 * b1 * b2 - a0 * b2 * b2 + a1 * b2 * b2 - a2 * b2 * b2) %
       iNormb) == 0);
}

bool RCyclotomic7Integer::isDivisbleBy(const AlgebraicInteger *ai) const {
  const RCyclotomic7Integer *rci(dynamic_cast<const RCyclotomic7Integer *>(ai));

  mpz_class a0(iC[0]), a1(iC[1]), a2(iC[2]);
  mpz_class b0(rci->iC[0]), b1(rci->iC[1]), b2(rci->iC[2]);

  mpz_class iNormb(rci->iNorm());

  return (
      ((-a0 * b0 * b0 - a1 * b0 * b0 + a2 * b0 * b0 - 2 * a0 * b0 * b1 +
        3 * a1 * b0 * b1 - a2 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 +
        3 * a0 * b0 * b2 - a2 * b0 * b2 + 2 * a0 * b1 * b2 - 3 * a1 * b1 * b2 +
        2 * a2 * b1 * b2 - 2 * a0 * b2 * b2 + 2 * a1 * b2 * b2 - a2 * b2 * b2) %
       iNormb) == 0 &&
      ((-a0 * b0 * b0 - 2 * a1 * b0 * b0 + 2 * a2 * b0 * b0 - a0 * b0 * b1 +
        3 * a1 * b0 * b1 + a0 * b1 * b1 - a1 * b1 * b1 - a2 * b1 * b1 +
        2 * a0 * b0 * b2 + 2 * a1 * b0 * b2 - 3 * a2 * b0 * b2 - a0 * b1 * b2 -
        2 * a1 * b1 * b2 + 3 * a2 * b1 * b2 + a1 * b2 * b2 - a2 * b2 * b2) %
       iNormb) == 0 &&
      ((-a0 * b0 * b0 + a2 * b0 * b0 - 3 * a0 * b0 * b1 + 2 * a1 * b0 * b1 +
        2 * a2 * b0 * b1 + 2 * a0 * b1 * b1 - a1 * b1 * b1 - 2 * a2 * b1 * b1 +
        3 * a0 * b0 * b2 - a1 * b0 * b2 - 2 * a2 * b0 * b2 - a1 * b1 * b2 +
        3 * a2 * b1 * b2 - a0 * b2 * b2 + a1 * b2 * b2 - a2 * b2 * b2) %
       iNormb) == 0);
}

void RCyclotomic7Integer::gcd(const AlgebraicInteger *ai) {
  if (ai->isEqualTo(0))
    return;

  const RCyclotomic7Integer *rci(dynamic_cast<const RCyclotomic7Integer *>(ai));

  if (this->isEqualTo(0)) {
    iC[0] = rci->iC[0];
    iC[1] = rci->iC[1];
    iC[2] = rci->iC[2];
    return;
  }

  mpz_class iGCD, iNorm1(abs(iNorm())), iNorm2(abs(rci->iNorm()));

  // -----------------------------------------
  // If one is isInvertible
  if (iNorm1 == 1 || iNorm2 == 1) {
    iC[0] = -1;
    iC[1] = -1;
    iC[2] = -1;
    return;
  }

  // -----------------------------------------
  // Some other trivial cases
  if (iC[0] == rci->iC[0] && iC[1] == rci->iC[1] && iC[2] == rci->iC[2])
    return;
  else if (rci->isDivisbleBy(this))
    return;
  else if (this->isDivisbleBy(rci)) {
    iC[0] = rci->iC[0];
    iC[1] = rci->iC[1];
    iC[2] = rci->iC[2];
    return;
  }

  // -----------------------------------------
  // GCD computations
  mpz_gcd(iGCD.get_mpz_t(), iNorm().get_mpz_t(), rci->iNorm().get_mpz_t());

  if (iGCD == 1) {
    iC[0] = -1;
    iC[1] = -1;
    iC[2] = -1;
    return;
  }

  // ---------------------------------------
  // We actually have to compute the GCD
  if (!iGCD.fits_slong_p())
    throw(string("RCyclotomic7Integer: gcd: Component too big"));

  auto iNormPrimesFactors(primeFactors<unsigned long>(abs(iGCD.get_si())));

  RCyclotomic7Integer rci1(*this), rci2(*rci);
  iC[0] = -1;
  iC[1] = -1;
  iC[2] = -1;

  for (auto iNPF : iNormPrimesFactors) {
    if (iNPF > iPrimesDecomposition_max)
      throw(string("RCyclotomic7Integer: gcd: Max prime number"));

    auto rciTemp(iPrimesDecomposition[iNPF]);

    for (auto rciFactor : rciTemp) {
      while (rci1.divideByIfDivisible(&rciFactor) &&
             rci2.divideByIfDivisible(&rciFactor))
        this->multiplyBy(&rciFactor);
    }
  }
}

bool RCyclotomic7Integer::isSquareOfIvertible() const {
  if (!isGreaterOEThan(0))
    return false;

  mpz_class uNorm(this->iNorm());
  if (uNorm != 1 && uNorm != -1)
    return false;

  return true;

  // TODO
}

bool RCyclotomic7Integer::isInvertible() const {
  return (abs(this->iNorm()) == 1);
}

vector<RCyclotomic7Integer> RCyclotomic7Integer::rciPrimeFactors() const {
  if ((iC[0] == 0 && iC[1] == 0 && iC[2] == 0) || isInvertible())
    return vector<RCyclotomic7Integer>(0);

  vector<RCyclotomic7Integer> rciPrimesFactors;

  mpz_class mpz_Norm(iNorm());
  if (!mpz_Norm.fits_slong_p())
    throw(string("RCyclotomic7Integer: rciPrimeFactors: Component too big"));

  auto iNormPrimesFactors(primeFactors<unsigned long>(abs(mpz_Norm.get_si())));

  for (auto iNPF : iNormPrimesFactors) {
    if (iNPF > iPrimesDecomposition_max)
      throw(string("RCyclotomic7Integer: rciPrimeFactors: Max prime number"));

    auto rciTemp(iPrimesDecomposition[iNPF]);

    for (auto rciF : rciTemp) {
      if (isDivisbleBy(
              &rciF)) // rciF divides either *this or the conjugate of *this
        rciPrimesFactors.push_back(rciF);
    }
  }

  return rciPrimesFactors;
}

map<RCyclotomic7Integer, unsigned int>
RCyclotomic7Integer::rciPrimeDecomposition() const {
  if ((iC[0] == 0 && iC[1] == 0 && iC[2] == 0) || isInvertible())
    return map<RCyclotomic7Integer, unsigned int>();

  map<RCyclotomic7Integer, unsigned int> rciDecomp;
  unsigned int iPower;

  mpz_class mpz_Norm(iNorm());
  if (!mpz_Norm.fits_slong_p())
    throw(string(
        "RCyclotomic7Integer: rciPrimeDecomposition: Component too big"));

  auto iNormPrimesFactors(primeFactors<unsigned long>(abs(mpz_Norm.get_si())));
  RCyclotomic7Integer rciTemp(*this);

  for (auto iNPF : iNormPrimesFactors) {
    if (iNPF > iPrimesDecomposition_max)
      throw(string("RCyclotomic7Integer: rciPrimeFactors: Max prime number"));

    auto rciFactors(iPrimesDecomposition[iNPF]);

    for (auto rciF : rciFactors) {
      iPower = 0;
      while (rciTemp.divideByIfDivisible(
          &rciF)) // rciF divides either *this or the conjugate of *this
        iPower++;

      if (iPower)
        rciDecomp[rciF] = iPower;
    }
  }

  return rciDecomp;
}

RCyclotomic7Integer &RCyclotomic7Integer::
operator=(const RCyclotomic7Integer &rci) {
  iC[0] = rci.iC[0];
  iC[1] = rci.iC[1];
  iC[2] = rci.iC[2];

  return *this;
}

RCyclotomic7Integer &RCyclotomic7Integer::
operator*=(const RCyclotomic7Integer &rci) {
  this->multiplyBy(&rci);
  return *this;
}

RCyclotomic7Integer &RCyclotomic7Integer::
operator/=(const RCyclotomic7Integer &rci) {
  this->divideBy(&rci);
  return *this;
}

RCyclotomic7Integer RCyclotomic7Integer::
operator*(const RCyclotomic7Integer &rci) const {
  RCyclotomic7Integer res(*this);
  res.multiplyBy(&rci);

  return res;
}

RCyclotomic7Integer RCyclotomic7Integer::
operator+(const RCyclotomic7Integer &rci) const {
  return RCyclotomic7Integer(
      {iC[0] + rci.iC[0], iC[1] + rci.iC[1], iC[2] + rci.iC[2]});
}

RCyclotomic7Integer RCyclotomic7Integer::
operator-(const RCyclotomic7Integer &rci) const {
  return RCyclotomic7Integer(
      {iC[0] - rci.iC[0], iC[1] - rci.iC[1], iC[2] - rci.iC[2]});
}

RCyclotomic7Integer RCyclotomic7Integer::operator-() const {
  return RCyclotomic7Integer({-iC[0], -iC[1], -iC[2]});
}

bool RCyclotomic7Integer::operator>(const RCyclotomic7Integer &rci) const {
  return !isLessOEThan(rci);
}

bool RCyclotomic7Integer::operator==(const long int &n) const {
  return (iC[0] == -n && iC[1] == -n && iC[2] == -n);
}

bool RCyclotomic7Integer::operator==(const RCyclotomic7Integer &rci) const {
  return (iC == rci.iC);
}

void RCyclotomic7Integer::removeSquareFactors() {
  map<RCyclotomic7Integer, unsigned int> rciDecomp(
      this->rciPrimeDecomposition());
  unsigned int i;
  for (auto it : rciDecomp) {
    unsigned int iMax(it.second % 2 ? it.second - 1 : it.second);
    for (i = 1; i <= iMax; i++)
      this->divideBy(&it.first);
  }
}

void RCyclotomic7Integer::set(AlgebraicInteger *ai) {
  iC = dynamic_cast<RCyclotomic7Integer *>(ai)->iC;
}

void RCyclotomic7Integer::set(const int &n) {
  iC = array<mpz_class, 3>({-n, -n, -n});
}

AlgebraicInteger *RCyclotomic7Integer::copy() const {
  return new RCyclotomic7Integer(*this);
}

AlgebraicInteger *RCyclotomic7Integer::aiCopyToInteger(const int &n) const {
  return new RCyclotomic7Integer(n);
}

std::string RCyclotomic7Integer::get_classname() const {
  return "RCyclotomic7Integer";
}
