#include "alvinfractions.h"

AlVinFractions::AlVinFractions(vector<AlgebraicInteger *> aiPossibleNorms2)
    : aiPossibleNorms2(aiPossibleNorms2),
      aiPossibleNorms2_max(aiPossibleNorms2[aiPossibleNorms2.size() - 1]),
      iLastMaximum(0), iBatchSize(1) {}

AlVinFractions::~AlVinFractions() {
  for (auto i : alvinfractions)
    delete i;

  for (auto i : aiPossibleNorms2)
    delete i;
}

vector<AlVinFraction *> AlVinFractions::getNextAlVinFraction() {
  if (!alvinfractions.size()) // Frist call of the function?
  {
    while (!alvinfractions.size())
      computeNextAlVinFractions(); // Compute next fractions
  }

  if (alvinfractions.size() &&
      alvinfractions_it ==
          alvinfractions.end()) // All the fractions were tested
  {
    for (auto i : alvinfractions)
      delete i;
    alvinfractions.clear();

    while (!alvinfractions.size())
      computeNextAlVinFractions(); // Compute next fractions
  }

  vector<AlVinFraction *> vf;

  AlgebraicInteger *aiNumerator((*alvinfractions_it)->aiNumerator);

  do // We may have to return more than one fraction
  {
    vf.push_back(*alvinfractions_it);

    ++alvinfractions_it;

  } while (alvinfractions_it != alvinfractions.end() &&
           aiNumerator->isEqualTo(*(*alvinfractions_it)->aiNumerator));

  return vf;
}

std::vector<AlgebraicInteger *> AlVinFractions::get_aiPossibleNorm2() const {
  return aiPossibleNorms2;
}

const std::vector<AlgebraicInteger *> *
AlVinFractions::get_ptraiPossibleNorm2() const {
  return &aiPossibleNorms2;
}
