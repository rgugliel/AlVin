#include "rationalinteger.h"

RationalInteger::~RationalInteger() {}

RationalInteger::RationalInteger() : val(0) {}

RationalInteger::RationalInteger(long int iVal) : val(iVal) {}

RationalInteger::RationalInteger(const RationalInteger &ri) : val(ri.val) {}

AlgebraicInteger *RationalInteger::copyToInteger(const int &n) const {
  return new RationalInteger(n);
}

AlgebraicInteger *RationalInteger::copy() const {
  return new RationalInteger(*this);
}

void RationalInteger::set(const int &n) { val = n; }

void RationalInteger::set(AlgebraicInteger *ai) {
  RationalInteger *ri(dynamic_cast<RationalInteger *>(ai));
  val = ri->val;
}

bool RationalInteger::isInvertible() const { return (val == 1 || val == -1); }

bool RationalInteger::isSquareOfIvertible() const { return (val == 1); }

void RationalInteger::gcd(const AlgebraicInteger *ai) {
  const RationalInteger *ri(dynamic_cast<const RationalInteger *>(ai));

  long u(val < 0 ? -val : val), v(ri->val < 0 ? -ri->val : ri->val);

  while (v != 0) {
    unsigned int r = u % v;
    u = v;
    v = r;
  }

  val = u;
}

long int RationalInteger::get_iValue() const { return val; }

void RationalInteger::removeSquareFactors() {
  unsigned long int iDivisor(3), iDivisorSquared;
  unsigned long int iResult(1);

  if (val == 1)
    return;

  while (!(val % 4))
    val /= 4;

  if (!(val % 2)) {
    iResult *= 2;
    do {
      val /= 2;
    } while (!(val % 2));
  }

  while (val > 1) {
    iDivisorSquared = iDivisor * iDivisor;
    while (!(val % iDivisorSquared))
      val /= iDivisorSquared;

    if (!(val % iDivisor)) {
      iResult *= iDivisor;
      do {
        val /= iDivisor;
      } while (!(val % iDivisor));
    }

    iDivisor += 2;
  }

  val = iResult;
}

bool operator<(const RationalInteger &i1, const RationalInteger &i2) {
  return (i1.val < i2.val);
}

bool RationalInteger::isLessThan(const AlgebraicInteger &ai) const {
  return (val < dynamic_cast<const RationalInteger &>(ai).val);
}

bool RationalInteger::isLessOEThan(const AlgebraicInteger &ai) const {
  return (val <= dynamic_cast<const RationalInteger &>(ai).val);
}

bool RationalInteger::isLessThan(const int &n) const { return (val < n); }

bool RationalInteger::isGreaterThan(const int &n) const { return (val > n); }

bool RationalInteger::isGreaterOEThan(const int &n) const { return (val >= n); }

bool RationalInteger::isEqualTo(const AlgebraicInteger &ai) const {
  return (val == dynamic_cast<const RationalInteger &>(ai).val);
}

bool RationalInteger::isEqualTo(const int &n) const { return (val == n); }

bool RationalInteger::isDivisibleBy(const AlgebraicInteger *ai) const {
  return !(val % dynamic_cast<const RationalInteger *>(ai)->val);
}

bool RationalInteger::divideByIfDivisible(const AlgebraicInteger *ai) {
  long int iV(dynamic_cast<const RationalInteger *>(ai)->val);

  if (val % iV)
    return false;

  val /= iV;

  return true;
}

void RationalInteger::divideBy(const AlgebraicInteger *ai) {
  val /= dynamic_cast<const RationalInteger *>(ai)->val;
}

void RationalInteger::add(const AlgebraicInteger *ai) {
  val += dynamic_cast<const RationalInteger *>(ai)->val;
}

void RationalInteger::multiplyBy(const int &n) { val *= n; }

void RationalInteger::multiplyBy(const AlgebraicInteger *ai) {
  val *= dynamic_cast<const RationalInteger *>(ai)->val;
}

void RationalInteger::substract(const AlgebraicInteger *ai) {
  val -= dynamic_cast<const RationalInteger *>(ai)->val;
}

void RationalInteger::opp() { val *= -1; }

ostream &RationalInteger::print(ostream &o) const {
  o << val;

  return o;
}

string RationalInteger::to_string(const string &strFormat,
                                  const bool &bProtect) const {
  return std::to_string(val);
}

double RationalInteger::to_double() const { return (double)val; }

string RationalInteger::get_classname() const { return "RationalInteger"; }

// ------------------------------------------------------------------------
// Operators
RationalInteger &RationalInteger::operator=(const RationalInteger &ri) {
  val = ri.val;

  return *this;
}

RationalInteger &RationalInteger::operator/=(const RationalInteger &ri) {
  val /= ri.val;
  return *this;
}

RationalInteger &RationalInteger::operator*=(const RationalInteger &ri) {
  val *= ri.val;
  return *this;
}

RationalInteger RationalInteger::operator-() const {
  return RationalInteger(-val);
}

RationalInteger RationalInteger::operator+(const RationalInteger &ri2) const {
  return RationalInteger(val + ri2.val);
}

RationalInteger RationalInteger::operator*(const RationalInteger &ri2) const {
  return RationalInteger(val * ri2.val);
}

RationalInteger RationalInteger::operator-(const RationalInteger &ri2) const {
  return RationalInteger(val - ri2.val);
}

bool RationalInteger::operator==(const long int &i) const { return (val == i); }

bool RationalInteger::operator==(const RationalInteger &ri) const {
  return (val == ri.val);
}

bool RationalInteger::operator>(const RationalInteger &ri) const {
  return !isLessOEThan(ri);
}

bool RationalInteger::operator<(const RationalInteger &ri) const {
  return isLessThan(ri);
}
