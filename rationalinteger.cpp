#include "rationalinteger.h"

RationalInteger::~RationalInteger() {}

RationalInteger::RationalInteger() : iVal(0) {}

RationalInteger::RationalInteger(long int iVal) : iVal(iVal) {}

RationalInteger::RationalInteger(const RationalInteger &ri) : iVal(ri.iVal) {}

AlgebraicInteger *RationalInteger::aiCopyToInteger(const int &n) const {
  return new RationalInteger(n);
}

AlgebraicInteger *RationalInteger::copy() const {
  return new RationalInteger(*this);
}

void RationalInteger::set(const int &n) { iVal = n; }

void RationalInteger::set(AlgebraicInteger *ai) {
  RationalInteger *ri(dynamic_cast<RationalInteger *>(ai));
  iVal = ri->iVal;
}

bool RationalInteger::isInvertible() const {
  return (iVal == 1 || iVal == -1);
}

bool RationalInteger::isSquareOfIvertible() const { return (iVal == 1); }

void RationalInteger::gcd(const AlgebraicInteger *ai) {
  const RationalInteger *ri(dynamic_cast<const RationalInteger *>(ai));

  long u(iVal < 0 ? -iVal : iVal), v(ri->iVal < 0 ? -ri->iVal : ri->iVal);

  while (v != 0) {
    unsigned int r = u % v;
    u = v;
    v = r;
  }

  iVal = u;
}

long int RationalInteger::get_iValue() const { return iVal; }

void RationalInteger::removeSquareFactors() {
  unsigned long int iDivisor(3), iDivisorSquared;
  unsigned long int iResult(1);

  if (iVal == 1)
    return;

  while (!(iVal % 4))
    iVal /= 4;

  if (!(iVal % 2)) {
    iResult *= 2;
    do {
      iVal /= 2;
    } while (!(iVal % 2));
  }

  while (iVal > 1) {
    iDivisorSquared = iDivisor * iDivisor;
    while (!(iVal % iDivisorSquared))
      iVal /= iDivisorSquared;

    if (!(iVal % iDivisor)) {
      iResult *= iDivisor;
      do {
        iVal /= iDivisor;
      } while (!(iVal % iDivisor));
    }

    iDivisor += 2;
  }

  iVal = iResult;
}

bool operator<(const RationalInteger &i1, const RationalInteger &i2) {
  return (i1.iVal < i2.iVal);
}

bool RationalInteger::isLessThan(const AlgebraicInteger &ai) const {
  return (iVal < dynamic_cast<const RationalInteger &>(ai).iVal);
}

bool RationalInteger::isLessOEThan(const AlgebraicInteger &ai) const {
  return (iVal <= dynamic_cast<const RationalInteger &>(ai).iVal);
}

bool RationalInteger::isLessThan(const int &n) const { return (iVal < n); }

bool RationalInteger::isGreaterThan(const int &n) const { return (iVal > n); }

bool RationalInteger::isGreaterOEThan(const int &n) const {
  return (iVal >= n);
}

bool RationalInteger::isEqualTo(const AlgebraicInteger &ai) const {
  return (iVal == dynamic_cast<const RationalInteger &>(ai).iVal);
}

bool RationalInteger::isEqualTo(const int &n) const { return (iVal == n); }

bool RationalInteger::isDivisbleBy(const AlgebraicInteger *ai) const {
  return !(iVal % dynamic_cast<const RationalInteger *>(ai)->iVal);
}

bool RationalInteger::divideByIfDivisible(const AlgebraicInteger *ai) {
  long int iV(dynamic_cast<const RationalInteger *>(ai)->iVal);

  if (iVal % iV)
    return false;

  iVal /= iV;

  return true;
}

void RationalInteger::divideBy(const AlgebraicInteger *ai) {
  iVal /= dynamic_cast<const RationalInteger *>(ai)->iVal;
}

void RationalInteger::add(const AlgebraicInteger *ai) {
  iVal += dynamic_cast<const RationalInteger *>(ai)->iVal;
}

void RationalInteger::multiplyBy(const int &n) { iVal *= n; }

void RationalInteger::multiplyBy(const AlgebraicInteger *ai) {
  iVal *= dynamic_cast<const RationalInteger *>(ai)->iVal;
}

void RationalInteger::substract(const AlgebraicInteger *ai) {
  iVal -= dynamic_cast<const RationalInteger *>(ai)->iVal;
}

void RationalInteger::opp() { iVal *= -1; }

ostream &RationalInteger::print(ostream &o) const {
  o << iVal;

  return o;
}

string RationalInteger::to_string(const string &strFormat,
                                  const bool &bProtect) const {
  return std::to_string(iVal);
}

double RationalInteger::to_double() const { return (double)iVal; }

string RationalInteger::get_classname() const { return "RationalInteger"; }

// ------------------------------------------------------------------------
// Operators
RationalInteger &RationalInteger::operator=(const RationalInteger &ri) {
  iVal = ri.iVal;

  return *this;
}

RationalInteger &RationalInteger::operator/=(const RationalInteger &ri) {
  iVal /= ri.iVal;
  return *this;
}

RationalInteger &RationalInteger::operator*=(const RationalInteger &ri) {
  iVal *= ri.iVal;
  return *this;
}

RationalInteger RationalInteger::operator-() const {
  return RationalInteger(-iVal);
}

RationalInteger RationalInteger::operator+(const RationalInteger &ri2) const {
  return RationalInteger(iVal + ri2.iVal);
}

RationalInteger RationalInteger::operator*(const RationalInteger &ri2) const {
  return RationalInteger(iVal * ri2.iVal);
}

RationalInteger RationalInteger::operator-(const RationalInteger &ri2) const {
  return RationalInteger(iVal - ri2.iVal);
}

bool RationalInteger::operator==(const long int &i) const {
  return (iVal == i);
}

bool RationalInteger::operator==(const RationalInteger &ri) const {
  return (iVal == ri.iVal);
}

bool RationalInteger::operator>(const RationalInteger &ri) const {
  return !isLessOEThan(ri);
}

bool RationalInteger::operator<(const RationalInteger &ri) const {
  return isLessThan(ri);
}
