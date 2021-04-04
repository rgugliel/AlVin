#include "rationalinteger_notreflective.h"

RationalInteger_NotReflective::RationalInteger_NotReflective(AlVin *v)
    : NotReflective(v) {}

void RationalInteger_NotReflective::createSystemEquations(
    NotReflective_Graph nrg) {
  cout << "\nSystem:" << endl;
  unsigned int iVariablesCount(nrg.iVariablesToCoeff.size());

  string strGlobalEquation("");

  for (auto aiNorm : aiPossibleNorm2) {
    // The cristallographic condition implies that the coefficient x_i is
    // divisible by (e,e) / gcd( (e,e), 2 * alpha_i ) Hence, we will do a
    // substitution.
    vector<AlgebraicInteger *> aiVariablesCoeff;
    for (unsigned int i(0); i < iVariablesCount; i++) {
      AlgebraicInteger *aiGCD(aiNorm->copy()), *aiRes(aiNorm->copy());
      aiGCD->gcd(ai2QF[nrg.iVariablesToCoeff[i]]);
      aiRes->divideBy(aiGCD);
      aiVariablesCoeff.push_back(aiRes);

      delete aiGCD;
    }

    string strEquation(""), strInequalities("x1 > 0");

    // ------------------------------------------------------------
    // norm equation
    for (unsigned int i(0); i < iVariablesCount; i++) {
      AlgebraicInteger *aiTemp(aiVariablesCoeff[i]->copy());

      aiTemp->multiplyBy(aiTemp);
      aiTemp->multiplyBy(nrg.aiVariablesCount[i]);
      strEquation += (i ? " + " : "-") + aiTemp->to_string() + " * x" +
                     to_string(i + 1) + "^2";

      if (nrg.iVariablesGreaterThan[i])
        strInequalities += " && x" + to_string(i) + " >= x" +
                           to_string(nrg.iVariablesGreaterThan[i]);

      delete aiTemp;
    }
    strEquation += " == " + aiNorm->to_string();

    // ----------------------------------------
    // other equations
    vector<AlgebraicInteger *> aiNumberVariablesCount;
    for (unsigned int i(0); i < iVariablesCount; i++)
      aiNumberVariablesCount.push_back(aiQF[0]->copyToInteger(0));
    AlgebraicInteger *aiTemp(aiQF[0]->copyToInteger(0));

    for (vector<short unsigned int>::const_iterator it(
             upper_bound(nrg.iGraphVertices.begin(), nrg.iGraphVertices.end(),
                         iDimension - 1));
         it != nrg.iGraphVertices.end(); ++it) {
      for (unsigned int i(0); i <= iDimension; i++) {
        if (aiVectors[*it][i] && nrg.iVariablesName[i]) {
          aiTemp->set(aiQF[i]);
          aiTemp->multiplyBy(aiVectors[*it][i]);
          aiTemp->multiplyBy(aiVariablesCoeff[nrg.iVariablesName[i] - 1]);

          aiNumberVariablesCount[nrg.iVariablesName[i] - 1]->add(aiTemp);
        }
      }

      strEquation += " && ";
      for (unsigned int i(0); i < iVariablesCount; i++)
        strEquation += (i ? " + " : "-") +
                       aiNumberVariablesCount[i]->to_string(strOFormat, true) +
                       "*x" + to_string(i + 1);
      strEquation += " == 0";

      for (unsigned int i(0); i < iVariablesCount; i++)
        aiNumberVariablesCount[i]->set(0);
    }

    // ----------------------------------------
    // inequalities
    strEquation = strInequalities + " && " + strEquation;

    strGlobalEquation += (strGlobalEquation == "" ? "" : " || ") + string("(") +
                         strEquation + ")";

    delete aiTemp;
    for (unsigned int j(0); j < iVariablesCount; j++) {
      delete aiVariablesCoeff[j];
      delete aiNumberVariablesCount[j];
    }
  }

  cout << "Reduce[{" << strGlobalEquation << "},{";
  for (unsigned int i(0); i < iVariablesCount; i++)
    cout << (i ? "," : "") << "x" << i + 1;
  cout << "},Integers]" << endl;
}
