#include "notreflective.h"

NotReflective::NotReflective(AlVin *v)
    : alvin(v), iDimension(v->get_iDimension()), aiQF(v->get_aiQF()),
      aiVectors(v->get_aiVectors()), strOFormat("mathematica"),
      strAlgebraicIntegerType(v->get_strAlgebraicIntegerType()),
      aiPossibleNorm2(v->get_aiPossibleNorm2()) {
  AlgebraicInteger *ai2(aiQF[0]->aiCopyToInteger(2));
  for (unsigned int i(0); i <= iDimension; i++) {
    ai2QF.push_back(aiQF[i]->copy());
    ai2QF[i]->multiplyBy(ai2);
  }
}

NotReflective::~NotReflective() {
  for (auto grs : graphs) {
    for (auto gr : grs) {
      for (auto c : gr.aiVariablesCount)
        delete c;
    }
  }
}

void NotReflective::Run() {
  prepareGraphsList();

  if (graphs.size())
    createSystemsEquations();
}

void NotReflective::prepareGraphsList() {
  // We get the euclidean graphs which cannot be extended
  CoxIter ci(alvin->get_iCoxeterMatrix(), iDimension);
  vector<vector<short unsigned int>> iGraphsNotExtendable(
      ci.bCanBeFiniteCovolume_complete());

  // If there is nothing to do
  if (!iGraphsNotExtendable.size()) {
    cout << "Try another method" << endl;
    return;
  }

  graphs = vector<vector<NotReflective_Graph>>(iDimension + 1);

  for (auto &iGVertices : iGraphsNotExtendable) {
    sort(iGVertices.begin(), iGVertices.end());
    if (iGVertices[0] >= iDimension) // Bad luck
      continue;

    vector<short unsigned int> iVariablesName;
    for (unsigned int i(0); i <= iDimension; i++)
      iVariablesName.push_back(i + 1);

    for (vector<short unsigned int>::reverse_iterator it(iGVertices.rbegin());
         it != iGVertices.rend(); ++it) {
      if (*it >= iDimension) // Intersting vertices are 0, ..., iDemension - 1
        continue;

      if (*it == iDimension - 1 || *aiQF[*it + 1] != *aiQF[*it + 2])
        iVariablesName[*it + 1] =
            0; // This coefficient of the vector is really 0
      else
        iVariablesName[*it + 1] = iVariablesName[*it + 2];
    }

    // ------------------------------------------------------------------------
    // We count the number of needed variable (we want to minimize this)
    vector<short unsigned int> iVariablesNameDistinct(iVariablesName);
    sort(iVariablesNameDistinct.begin(), iVariablesNameDistinct.end());
    iVariablesNameDistinct = vector<short unsigned int>(
        iVariablesNameDistinct.begin(),
        unique(iVariablesNameDistinct.begin(), iVariablesNameDistinct.end()));
    unsigned int iVariablesCount(iVariablesNameDistinct[0] == 0
                                     ? iVariablesNameDistinct.size() - 1
                                     : iVariablesNameDistinct.size());

    // ------------------------------------------------------------------------
    // Renaming of the variables, order of the variable and counting
    unsigned int iLastVariable(0), iVariableIndex(0);
    vector<short unsigned int> iVariablesGreaterThan(iVariablesCount, 0);
    vector<AlgebraicInteger *> aiVariableCount;
    vector<short unsigned int> iVariablesToCoeff(iVariablesCount, 0);

    for (unsigned int i(0); i < iVariablesCount; i++)
      aiVariableCount.push_back(aiQF[0]->aiCopyToInteger(0));

    for (unsigned int i(0); i <= iDimension; i++) {
      if (iVariablesName[i]) {
        if (iVariablesName[i] != iLastVariable) {
          if (iVariableIndex && aiQF[i]->bIsEqualTo(*aiQF[i - 1]))
            iVariablesGreaterThan[iVariableIndex] = iVariableIndex + 1;

          iVariablesToCoeff[iVariableIndex] = i;
          iVariableIndex++;
        }

        iLastVariable = iVariablesName[i];
        iVariablesName[i] = iVariableIndex;

        aiVariableCount[iVariableIndex - 1]->add(aiQF[i]);
      }
    }

    graphs[iVariablesCount].push_back(
        NotReflective_Graph({iGVertices, iVariablesName, iVariablesToCoeff,
                             iVariablesGreaterThan, aiVariableCount}));
  }
}

void NotReflective::createSystemsEquations() {
  for (unsigned int i(2); i <= (5 <= iDimension ? 5 : iDimension); i++) {
    for (auto nrg : graphs[i])
      createSystemEquations(nrg);
  }
}
