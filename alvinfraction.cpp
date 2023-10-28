#include "alvinfraction.h"

AlVinFraction::AlVinFraction(AlgebraicInteger *aiX0, AlgebraicInteger *aiNorm2,
                             AlgebraicInteger *aiNumerator)
    : aiX0(aiX0), aiNorm2(aiNorm2), aiNumerator(aiNumerator) {}

AlVinFraction::~AlVinFraction() {
  delete aiX0;
  delete aiNorm2;
  delete aiNumerator;
}

bool operator<(const AlVinFraction &f1, const AlVinFraction &f2) {
  return (f1.aiNumerator->isLessThan(*f2.aiNumerator));
}

bool AlVinFraction::operator==(const AlVinFraction &f2) const {
  return (aiX0->isEqualTo(*f2.aiX0) && aiNorm2->isEqualTo(*f2.aiNorm2));
}

bool AlVinFraction::operator!=(const AlVinFraction &f2) const {
  return (!aiX0->isEqualTo(*f2.aiX0) || !aiNorm2->isEqualTo(*f2.aiNorm2));
}

ostream &operator<<(ostream &o, AlVinFraction const &vf) {
  o << *vf.aiX0 << "^2"
    << " / " << *vf.aiNorm2;

  return o;
}

bool isLessThanPtrAlVinFraction(const AlVinFraction *a,
                                const AlVinFraction *b) {
  return (*a < *b);
}
