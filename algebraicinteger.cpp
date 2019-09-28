#include "algebraicinteger.h"

bool isLessThanPtrAlgebraicInteger(AlgebraicInteger *a, AlgebraicInteger *b) {
  return (a->bIsLessThan(*b));
}

bool isEqualToPtrAlgebraicInteger(AlgebraicInteger *a, AlgebraicInteger *b) {
  return (a->bIsEqualTo(*b));
}

AlgebraicInteger::~AlgebraicInteger() {}

ostream &AlgebraicInteger::print(ostream &o) const { return o; }

bool AlgebraicInteger::operator==(const AlgebraicInteger &ai) const {
  return bIsEqualTo(ai);
}

bool AlgebraicInteger::operator!=(const int &n) const { return !bIsEqualTo(n); }

bool AlgebraicInteger::operator!=(const AlgebraicInteger &ai) const {
  return !bIsEqualTo(ai);
}

ostream &operator<<(ostream &o, AlgebraicInteger const &i) {
  i.print(o);

  return o;
}