#include "rationalinteger_notreflective.h"

RationalInteger_NotReflective::RationalInteger_NotReflective(AlVin *v)
    : NotReflective(v) {}

void RationalInteger_NotReflective::createSystemEquations(
    NotReflective_Graph nrg) {
  cout << "\nSystem:" << endl;
  unsigned int iVariablesCount(nrg.iVariablesToCoeff.size());

  string strGlobalEquation("");

  for (auto norm2 : possibleNorm2) {
    // The cristallographic condition implies that the coefficient x_i is
    // divisible by (e,e) / gcd( (e,e), 2 * alpha_i ) Hence, we will do a
    // substitution.
    vector<AlgebraicInteger *> variablesCoeff;
    for (unsigned int i(0); i < iVariablesCount; i++) {
      AlgebraicInteger *gcd(norm2->copy()), *res(norm2->copy());
      gcd->gcd(ai2QF[nrg.iVariablesToCoeff[i]]);
      res->divideBy(gcd);
      variablesCoeff.push_back(res);

      delete gcd;
    }

    string strEquation(""), strInequalities("x1 > 0");

    // ------------------------------------------------------------
    // norm equation
    for (unsigned int i(0); i < iVariablesCount; i++) {
      AlgebraicInteger *temp(variablesCoeff[i]->copy());

      temp->multiplyBy(temp);
      temp->multiplyBy(nrg.variablesCount[i]);
      strEquation += (i ? " + " : "-") + temp->to_string() + " * x" +
                     to_string(i + 1) + "^2";

      if (nrg.variablesGreaterThan[i])
        strInequalities += " && x" + to_string(i) + " >= x" +
                           to_string(nrg.variablesGreaterThan[i]);

      delete temp;
    }
    strEquation += " == " + norm2->to_string();

    // ----------------------------------------
    // other equations
    vector<AlgebraicInteger *> numberVariablesCount;
    for (unsigned int i(0); i < iVariablesCount; i++)
      numberVariablesCount.push_back(qf[0]->copyToInteger(0));
    AlgebraicInteger *temp(qf[0]->copyToInteger(0));

    for (vector<short unsigned int>::const_iterator it(
             upper_bound(nrg.iGraphVertices.begin(), nrg.iGraphVertices.end(),
                         iDimension - 1));
         it != nrg.iGraphVertices.end(); ++it) {
      for (unsigned int i(0); i <= iDimension; i++) {
        if (vectors[*it][i] && nrg.iVariablesName[i]) {
          temp->set(qf[i]);
          temp->multiplyBy(vectors[*it][i]);
          temp->multiplyBy(variablesCoeff[nrg.iVariablesName[i] - 1]);

          numberVariablesCount[nrg.iVariablesName[i] - 1]->add(temp);
        }
      }

      strEquation += " && ";
      for (unsigned int i(0); i < iVariablesCount; i++)
        strEquation += (i ? " + " : "-") +
                       numberVariablesCount[i]->to_string(strOFormat, true) +
                       "*x" + to_string(i + 1);
      strEquation += " == 0";

      for (unsigned int i(0); i < iVariablesCount; i++)
        numberVariablesCount[i]->set(0);
    }

    // ----------------------------------------
    // inequalities
    strEquation = strInequalities + " && " + strEquation;

    strGlobalEquation += (strGlobalEquation == "" ? "" : " || ") + string("(") +
                         strEquation + ")";

    delete temp;
    for (unsigned int j(0); j < iVariablesCount; j++) {
      delete variablesCoeff[j];
      delete numberVariablesCount[j];
    }
  }

  cout << "Reduce[{" << strGlobalEquation << "},{";
  for (unsigned int i(0); i < iVariablesCount; i++)
    cout << (i ? "," : "") << "x" << i + 1;
  cout << "},Integers]" << endl;
}
