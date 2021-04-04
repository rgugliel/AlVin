#include "alvinfraction.h"

AlVinFraction::AlVinFraction(AlgebraicInteger *aiX0, AlgebraicInteger *norm2,
                             AlgebraicInteger *numerator)
    : aiX0(aiX0), norm2(norm2), numerator(numerator) {}

AlVinFraction::~AlVinFraction() {
  delete aiX0;
  delete norm2;
  delete numerator;
}

bool operator<(const AlVinFraction &f1, const AlVinFraction &f2) {
  return (f1.numerator->isLessThan(*f2.numerator));
}

bool AlVinFraction::operator==(const AlVinFraction &f2) const {
  return (aiX0->isEqualTo(*f2.aiX0) && norm2->isEqualTo(*f2.norm2));
}

bool AlVinFraction::operator!=(const AlVinFraction &f2) const {
  return (!aiX0->isEqualTo(*f2.aiX0) || !norm2->isEqualTo(*f2.norm2));
}

ostream &operator<<(ostream &o, AlVinFraction const &vf) {
  o << *vf.aiX0 << "^2"
    << " / " << *vf.norm2;

  return o;
}

bool isLessThanPtrAlVinFraction(const AlVinFraction *a,
                                const AlVinFraction *b) {
  return (*a < *b);
}
