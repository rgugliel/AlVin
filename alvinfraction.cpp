#include "alvinfraction.h"

AlVinFraction::AlVinFraction(AlgebraicInteger *x0, AlgebraicInteger *norm2,
                             AlgebraicInteger *numerator)
    : x0(x0), norm2(norm2), numerator(numerator) {}

AlVinFraction::~AlVinFraction() {
  delete x0;
  delete norm2;
  delete numerator;
}

bool operator<(const AlVinFraction &f1, const AlVinFraction &f2) {
  return (f1.numerator->isLessThan(*f2.numerator));
}

bool AlVinFraction::operator==(const AlVinFraction &f2) const {
  return (x0->isEqualTo(*f2.x0) && norm2->isEqualTo(*f2.norm2));
}

bool AlVinFraction::operator!=(const AlVinFraction &f2) const {
  return (!x0->isEqualTo(*f2.x0) || !norm2->isEqualTo(*f2.norm2));
}

ostream &operator<<(ostream &o, AlVinFraction const &vf) {
  o << *vf.x0 << "^2"
    << " / " << *vf.norm2;

  return o;
}

bool isLessThanPtrAlVinFraction(const AlVinFraction *a,
                                const AlVinFraction *b) {
  return (*a < *b);
}
