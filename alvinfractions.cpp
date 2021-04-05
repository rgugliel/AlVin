#include "alvinfractions.h"

AlVinFractions::AlVinFractions(vector<AlgebraicInteger *> possibleNorms2)
    : possibleNorms2(possibleNorms2),
      possibleNorms2Max(possibleNorms2[possibleNorms2.size() - 1]),
      lastMaximum(0), batchSize(1) {}

AlVinFractions::~AlVinFractions() {
  for (auto i : alvinFractions)
    delete i;

  for (auto i : possibleNorms2)
    delete i;
}

vector<AlVinFraction *> AlVinFractions::getNextAlVinFraction() {
  if (!alvinFractions.size()) // Frist call of the function?
  {
    while (!alvinFractions.size())
      computeNextAlVinFractions(); // Compute next fractions
  }

  if (alvinFractions.size() &&
      alvinfractions_it ==
          alvinFractions.end()) // All the fractions were tested
  {
    for (auto i : alvinFractions)
      delete i;
    alvinFractions.clear();

    while (!alvinFractions.size())
      computeNextAlVinFractions(); // Compute next fractions
  }

  vector<AlVinFraction *> vf;

  AlgebraicInteger *numerator((*alvinfractions_it)->numerator);

  do // We may have to return more than one fraction
  {
    vf.push_back(*alvinfractions_it);

    ++alvinfractions_it;

  } while (alvinfractions_it != alvinFractions.end() &&
           numerator->isEqualTo(*(*alvinfractions_it)->numerator));

  return vf;
}

std::vector<AlgebraicInteger *> AlVinFractions::get_possibleNorm2() const {
  return possibleNorms2;
}

const std::vector<AlgebraicInteger *> *
AlVinFractions::get_ptrPossibleNorm2() const {
  return &possibleNorms2;
}
