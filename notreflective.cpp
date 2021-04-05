#include "notreflective.h"

NotReflective::NotReflective(AlVin *v)
    : alvin(v), iDimension(v->get_dimension()), qf(v->get_qf()),
      vectors(v->get_vectors()), strOFormat("mathematica"),
      strAlgebraicIntegerType(v->get_strAlgebraicIntegerType()),
      possibleNorm2(v->get_possibleNorm2()) {
  AlgebraicInteger *ai2(qf[0]->copyToInteger(2));
  for (unsigned int i(0); i <= iDimension; i++) {
    ai2QF.push_back(qf[i]->copy());
    ai2QF[i]->multiplyBy(ai2);
  }
}

NotReflective::~NotReflective() {
  for (auto grs : graphs) {
    for (auto gr : grs) {
      for (auto c : gr.variablesCount)
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
  CoxIter ci(alvin->get_coxeterMatrix(), iDimension);
  vector<vector<short unsigned int>> iGraphsNotExtendable(
      ci.canBeFiniteCovolume_complete());

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

    vector<short unsigned int> variablesName;
    for (unsigned int i(0); i <= iDimension; i++)
      variablesName.push_back(i + 1);

    for (vector<short unsigned int>::reverse_iterator it(iGVertices.rbegin());
         it != iGVertices.rend(); ++it) {
      if (*it >= iDimension) // Intersting vertices are 0, ..., iDemension - 1
        continue;

      if (*it == iDimension - 1 || *qf[*it + 1] != *qf[*it + 2])
        variablesName[*it + 1] =
            0; // This coefficient of the vector is really 0
      else
        variablesName[*it + 1] = variablesName[*it + 2];
    }

    // ------------------------------------------------------------------------
    // We count the number of needed variable (we want to minimize this)
    vector<short unsigned int> variablesNameDistinct(variablesName);
    sort(variablesNameDistinct.begin(), variablesNameDistinct.end());
    variablesNameDistinct = vector<short unsigned int>(
        variablesNameDistinct.begin(),
        unique(variablesNameDistinct.begin(), variablesNameDistinct.end()));
    unsigned int variablesCount(variablesNameDistinct[0] == 0
                                    ? variablesNameDistinct.size() - 1
                                    : variablesNameDistinct.size());

    // ------------------------------------------------------------------------
    // Renaming of the variables, order of the variable and counting
    unsigned int iLastVariable(0), iVariableIndex(0);
    vector<short unsigned int> iVariablesGreaterThan(variablesCount, 0);
    vector<AlgebraicInteger *> variableCount;
    vector<short unsigned int> iVariablesToCoeff(variablesCount, 0);

    for (unsigned int i(0); i < variablesCount; i++)
      variableCount.push_back(qf[0]->copyToInteger(0));

    for (unsigned int i(0); i <= iDimension; i++) {
      if (variablesName[i]) {
        if (variablesName[i] != iLastVariable) {
          if (iVariableIndex && qf[i]->isEqualTo(*qf[i - 1]))
            iVariablesGreaterThan[iVariableIndex] = iVariableIndex + 1;

          iVariablesToCoeff[iVariableIndex] = i;
          iVariableIndex++;
        }

        iLastVariable = variablesName[i];
        variablesName[i] = iVariableIndex;

        variableCount[iVariableIndex - 1]->add(qf[i]);
      }
    }

    graphs[variablesCount].push_back(
        NotReflective_Graph({iGVertices, variablesName, iVariablesToCoeff,
                             iVariablesGreaterThan, variableCount}));
  }
}

void NotReflective::createSystemsEquations() {
  for (unsigned int i(2); i <= (5 <= iDimension ? 5 : iDimension); i++) {
    for (auto nrg : graphs[i])
      createSystemEquations(nrg);
  }
}
