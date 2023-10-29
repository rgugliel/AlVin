#include "algebraicinteger.h"

bool isLessThanPtrAlgebraicInteger(AlgebraicInteger *a, AlgebraicInteger *b) {
  return (a->isLessThan(*b));
}

bool isEqualToPtrAlgebraicInteger(AlgebraicInteger *a, AlgebraicInteger *b) {
  return (a->isEqualTo(*b));
}

AlgebraicInteger::~AlgebraicInteger() {}

ostream &AlgebraicInteger::print(ostream &o) const { return o; }

bool AlgebraicInteger::operator==(const AlgebraicInteger &ai) const {
  return isEqualTo(ai);
}

bool AlgebraicInteger::operator!=(const int &n) const { return !isEqualTo(n); }

bool AlgebraicInteger::operator!=(const AlgebraicInteger &ai) const {
  return !isEqualTo(ai);
}

ostream &operator<<(ostream &o, AlgebraicInteger const &i) {
  i.print(o);

  return o;
}