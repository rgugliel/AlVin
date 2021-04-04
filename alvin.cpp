#include "alvin.h"

AlVin::AlVin(const string &strOuputMathematicalFormat, const bool &bWriteInfo,
             const bool &bDebug)
    : bComputeInvariantsPolyhedron(false), bWriteInfo(bWriteInfo),
      bDebug(bDebug), iCreateImage(-1), vectorsCountSecond(0), vf(nullptr),
      ptrCI(nullptr), strOuputMathematicalFormat(strOuputMathematicalFormat) {}

AlVin::~AlVin() {
  for (auto i : qf)
    delete i;

  for (auto i : bilinearProducts)
    delete i;

  for (auto v : vectors) {
    for (auto i : v)
      delete i;
  }

  if (vf != nullptr)
    delete vf;

  for (auto i : candidateVectors) {
    for (auto j : i)
      delete j;
  }

  if (ptrCI != nullptr)
    delete ptrCI;
}

void AlVin::initializations() {
  if (qf.size() < 3)
    throw(string("Quadratic form has not enough coefficients"));

  ptrCI = nullptr;
  bilinearProducts.clear();
  dimension = qf.size() - 1;

  for (auto i : qf)
    i->removeSquareFactors();

  sort(qf.begin() + 1, qf.end(), isLessThanPtrAlgebraicInteger);

  // ------------------------------------------------
  // Test: quotients of the positive elements of the QF should not be the square
  // of an invertible element
  for (unsigned int i(1); i <= dimension; i++) {
    while (i < dimension && *qf[i + 1] == *qf[i])
      i++;

    for (unsigned int j(i + 1); j <= dimension; j++) {
      if (qf[j]->isDivisibleBy(qf[i])) {
        unique_ptr<AlgebraicInteger> aiTest(qf[j]->copy());
        aiTest->divideBy(qf[i]);

        if (aiTest->isSquareOfIvertible() && get_strField() == "RC7") {
          cout << "Check that " << *aiTest
               << " is not the square of an invertible element" << endl;
          strFinalInformation += "Check that " + aiTest->to_string() +
                                 " is not the square of an invertible element";
        }
      }
    }
  }

  // ------------------------------------------------
  // Test: GCD
  unique_ptr<AlgebraicInteger> igcd(qf[0]->copy());
  for (unsigned int i(1); i <= dimension; i++)
    igcd->gcd(qf[i]);

  if (!igcd->isInvertible())
    throw(
        string("GCD of the coefficients of the quadratic form should be one"));

  // ------------------------------------------------
  // Other stuff
  strAlgebraicIntegerType = qf[0]->get_classname();

  findFirstVectors();
}

void AlVin::findFirstVectors() {
  vectorsCount = dimension;
  vectorsCountSecond = 0;

  // ------------------------------------------------
  // We count the multiplicity of the coefficients
  unsigned int iCount(1);

  unique_ptr<AlgebraicInteger> aiLastCoeff(qf[1]->copy());

  for (unsigned i(2); i <= dimension; i++) {
    if (*qf[i] != *aiLastCoeff) {
      iQBlocksSize.push_back(iCount);
      iCount = 0;
    }

    aiLastCoeff->set(qf[i]);
    iCount++;
  }

  iQBlocksSize.push_back(iCount);

  // ------------------------------------------------
  // Vectors and gram matrix
  coxeterMatrix = vector<vector<unsigned int>>(
      dimension, vector<unsigned int>(dimension, 2));
  componentLessThan = vector<unsigned int>(dimension + 1, 0);

  unsigned int iVectorIndex(1), i;
  for (auto iBlockSize : iQBlocksSize) {
    for (i = 1; i <= iBlockSize; i++) {
      vector<AlgebraicInteger *> aiV;
      for (unsigned int j(0); j <= dimension; j++)
        aiV.push_back(qf[0]->copyToInteger(0));

      if (i >= 2)
        componentLessThan[iVectorIndex] = iVectorIndex - 1;

      if (i == iBlockSize)
        aiV[iVectorIndex]->set(-1);
      else {
        aiV[iVectorIndex]->set(-1);
        aiV[iVectorIndex + 1]->set(1);

        coxeterMatrix[iVectorIndex - 1][iVectorIndex] =
            coxeterMatrix[iVectorIndex][iVectorIndex - 1] =
                i < iBlockSize - 1 ? 3 : 4;
      }

      vectors.push_back(aiV);
      addVectorChild(aiV);

      iVectorIndex++;
    }
  }
}

bool AlVin::Run(unsigned int iMinVectors, unsigned int iMaxVectors,
                bool bLastCheckFV) {
  if (!PreRun())
    return false;

  unsigned int iCandidatesCount, i, j;

  while ((!iMaxVectors || vectorsCount < iMaxVectors)) {
    vector<AlVinFraction *> vfs(vf->getNextAlVinFraction());
    for (auto frac : vfs)
      findVector(frac->aiX0, frac->norm2);

    // ------------------------------------------------
    // Are the candidates compatible between themselves
    iCandidatesCount = candidateVectors.size();

    for (i = 0; i < iCandidatesCount; i++) {
      for (j = i + 1; j < iCandidatesCount; j++) {
        unique_ptr<AlgebraicInteger> aiProd(
            bilinearProduct(candidateVectors[i], candidateVectors[j]));

        if (aiProd->isGreaterThan(0)) // This should not happen
          throw(string("AlVin::Run: Incompatible candidates"));
      }
    }

    // ------------------------------------------------
    // Adding the vectors
    for (i = 0;
         i < iCandidatesCount && (!iMaxVectors || vectorsCount < iMaxVectors);
         i++) {
      if (bWriteInfo) {
        printFoundVector(candidateVectors[i], vectorsCount + 1, false);
      }

      try {
        addVector(candidateVectors[i]);
        addVectorChild(candidateVectors[i]);
      } catch (string &strE) {
        candidateVectors[i].clear();
        throw(strE);
      }

      // So that we don't try do free the memory (will be done via aiVectors)
      candidateVectors[i].clear();

      // --------------------------------------------------
      // Test volume
      if (iMinVectors && vectorsCount < iMinVectors)
        continue; // No covolume test

      ptrCI = new CoxIter(coxeterMatrix, dimension);
      ptrCI->set_checkCofiniteness(true);

      if (!ptrCI->canBeFiniteCovolume()) {
        delete ptrCI;
        ptrCI = nullptr;
        continue;
      }

      if ((vectorsCount == iMaxVectors && !bLastCheckFV) ||
          ptrCI->checkCovolumeFiniteness() == 1) {
        if (bWriteInfo) {
          print_finallInformation();

          cout << "Algorithm ended\n" << endl;
        }

        string strFilename(to_string(dimension) + "-");
        for (j = 0; j <= dimension; j++)
          strFilename +=
              qf[j]->to_string("filename") + (j < dimension ? "," : "");

        if (!ptrCI->writeGraph("output/" + strFilename) ||
            !ptrCI->writeGraphToDraw("output/" + strFilename))
          cout << "Error:\n\tCheck that the folder 'output/' exists and is "
                  "writable"
               << endl;
        else if (bWriteInfo) {
          cout << "Graph written in the file: \n\t"
               << "output/" + strFilename << ".coxiter" << endl;

          if (bWriteInfo) {
            if (iCreateImage == 1 ||
                (iCreateImage == -1 && vectorsCount <= 25)) {
              FILE *fin;
              string strCmd("dot -Tjpg -ooutput/" + strFilename +
                            ".jpg  output/" + strFilename + ".graphviz");
              if ((fin = popen(strCmd.c_str(), "r")))
                pclose(fin);
            } else {
              cout << "\nCommand to create the image: \n\tdot -Tjpg -ooutput/"
                   << strFilename << ".jpg  output/" << strFilename
                   << ".graphviz\n"
                   << endl;
            }
          }
        }

        return true;
      }

      if (ptrCI != nullptr) {
        delete ptrCI;
        ptrCI = nullptr;
      }
    }

    candidateVectors.clear();
  }

  CoxIter ci(coxeterMatrix, dimension);

  string strFilename("nt-" + to_string(dimension) + "-");
  for (j = 0; j <= dimension; j++)
    strFilename += qf[j]->to_string() + (j < dimension ? "," : "");

  ci.writeGraph("output/" + strFilename);
  ci.writeGraphToDraw("output/" + strFilename);

  if (bWriteInfo) {
    if (iCreateImage == 1 || (iCreateImage == -1 && vectorsCount <= 25)) {
      FILE *fin;
      string strCmd("dot -Tjpg -ooutput/" + strFilename + ".jpg output/" +
                    strFilename + ".graphviz");
      if ((fin = popen(strCmd.c_str(), "r")))
        pclose(fin);
    } else
      cout << "\nCommand to create the image: \n\tdot -Tjpg -ooutput/"
           << strFilename << ".jpg output/" << strFilename << ".graphviz\n"
           << endl;

    print_finallInformation();
  }

  return false;
}

void AlVin::addVector(const vector<AlgebraicInteger *> &iVect) {
  int iWeight;

  AlgebraicInteger *aiNorm(bilinearProduct(iVect, iVect));

  coxeterMatrix.push_back(vector<unsigned int>(vectorsCount + 1, 2));

  vectors.push_back(iVect);
  for (unsigned int i(0); i < vectorsCount; i++) {
    AlgebraicInteger *aiProd(bilinearProduct(vectors[i], iVect));

    iWeight = -2;

    if (aiProd->isEqualTo(0))
      iWeight = 2;
    else {
      AlgebraicInteger *numerator(aiProd->copy());
      numerator->multiplyBy(aiProd);

      AlgebraicInteger *denominator(bilinearProduct(vectors[i], vectors[i]));
      denominator->multiplyBy(aiNorm);

      if (denominator->isLessThan(*numerator))
        iWeight = 1; // Dotted edge
      else {
        AlgebraicInteger *gcd(numerator->copy());
        gcd->gcd(denominator);

        numerator->divideBy(gcd);
        denominator->divideBy(gcd);
        delete gcd;

        if (numerator->isEqualTo(1)) {
          if (denominator->isEqualTo(4))
            iWeight = 3;
          else if (denominator->isEqualTo(2))
            iWeight = 4;
          else if (denominator->isEqualTo(1))
            iWeight = 0; // Infty
        } else if (numerator->isEqualTo(3)) {
          if (denominator->isEqualTo(4))
            iWeight = 6;
        }

        if (iWeight == -2)
          iWeight = addVector_findWeight(numerator, denominator);
      }

      if (iWeight == -2)
        throw(string("Unknown weight betweens vectors " + to_string(i + 1) +
                     " and " + to_string(vectorsCount + 1)));

      delete denominator;
      delete numerator;
    }

    delete aiProd;

    coxeterMatrix[i].push_back(iWeight);
    coxeterMatrix[vectorsCount][i] = iWeight;
  }

  vectorsCount++;
  vectorsCountSecond++;

  delete aiNorm;
}

void AlVin::printFoundVector(std::vector<AlgebraicInteger *> aiV,
                             const unsigned int &iIndex,
                             const bool &bFirst) const {
  if (strOuputMathematicalFormat == "latex") {
    cout << "\te_" << (iIndex > 9 ? "{" : "") << iIndex
         << (iIndex > 9 ? "}" : "") << " = \\left(";
    for (auto it(aiV.begin()); it != aiV.end(); ++it)
      cout << (it == aiV.begin() ? "" : ", ") << (*it)->to_string("latex");
    cout << "\\right)" << endl;
  } else if (strOuputMathematicalFormat == "mathematica") {
    if (bFirst) {
      cout << "vectors := {{";
      for (auto it(aiV.begin()); it != aiV.end(); ++it)
        cout << (it == aiV.begin() ? "" : ", ")
             << (*it)->to_string("mathematica");
      cout << "}};" << endl;
    } else {
      cout << "AppendTo[vectors, {";
      for (auto it(aiV.begin()); it != aiV.end(); ++it)
        cout << (it == aiV.begin() ? "" : ", ")
             << (*it)->to_string("mathematica");
      cout << "}]; (* vector " << iIndex << " *); " << endl;
    }
  } else if (strOuputMathematicalFormat == "pari") {
    if (bFirst)
      cout << "vectors = Mat( [";
    else
      cout << "vectors = concat(vectors, [";

    for (auto it(aiV.begin()); it != aiV.end(); ++it)
      cout << (it == aiV.begin() ? "" : ", ") << (*it)->to_string("pari");
    cout << "] );" << endl;
  } else {
    cout << "\te" << iIndex << " = (";
    for (auto it(aiV.begin()); it != aiV.end(); ++it)
      cout << (it == aiV.begin() ? "" : ", ") << **it;
    cout << ")" << endl;
  }
}

void AlVin::print_iQF() const {
  cout << "Quadratic form (" << dimension << ",1): ";
  for (auto it(qf.begin()); it != qf.end(); it++) {
    if (it == qf.begin()) {
      unique_ptr<AlgebraicInteger> ptr((*it)->copy());
      ptr->opp();
      cout << *ptr;
    } else
      cout << ", " << **it;
  }
  cout << endl;

  cout << "Field of definition: " << get_strField() << endl;
}

void AlVin::print_initialInformationChild() const {}

void AlVin::print_initialInformation() const {
  print_iQF();

  if (bDebug) {
    cout << "Conditions on the coefficients: \n";
    for (unsigned int i(1); i <= dimension; i++) {
      if (componentLessThan[i])
        cout << "\tx" << i << " <= x" << componentLessThan[i] << endl;
    }
    cout << endl;

    const vector<AlgebraicInteger *> *ptraiPN2(vf->get_ptraiPossibleNorm2());
    unsigned int iMax(ptraiPN2->size());

    cout << "Possible values for (e,e): ";
    for (unsigned int i(0); i < iMax; i++)
      cout << (i ? ", " : "") << *(*ptraiPN2)[i];
    cout << endl;
  }

  cout << endl;

  print_initialInformationChild();
  if (strOuputMathematicalFormat == "mathematica") // quadratic form
  {
    cout << "F[x_, y_] := Sum[x[[i]]*y[[i]]*qform[[i]], {i, 1, "
            "Length[x]}];\nS[x_, y_] := F[x, y]/Sqrt[F[x, x]*F[y, y]];"
         << endl;
    cout << "qform := {";
    for (unsigned int i(0); i <= dimension; i++)
      cout << (i ? ", " : "-(") << qf[i]->to_string("mathematica")
           << (i ? "" : ")");
    cout << "};" << endl;
  } else if (strOuputMathematicalFormat == "pari") {
    cout << "qform = matdiagonal([";
    for (unsigned int i(0); i <= dimension; i++)
      cout << (i ? ", " : "-(") << qf[i]->to_string("pari") << (i ? "" : ")");
    cout << "]);" << endl;

    cout << "S = (x,y) -> qfeval(qform,x,y)/sqrt( qfeval(qform,x) * "
            "qfeval(qform,y) );"
         << endl;
  } else
    cout << "Vectors: " << endl;

  for (unsigned int i(0); i < vectorsCount; i++)
    printFoundVector(vectors[i], i + 1, i == 0);
}

void AlVin::print_finallInformation() const {
  if (strOuputMathematicalFormat == "mathematica") // quadratic form
  {
    cout << "grammatrix := Table[S[vectors[[i]], vectors[[j]]], {i, 1, "
            "Length[vectors]}, {j, 1, Length[vectors]}];"
         << endl;
  } else if (strOuputMathematicalFormat == "pari") {
    cout << "grammatrix = matrix( " << vectorsCount << "," << vectorsCount
         << ", i, j, S(vectors[i,], vectors[j,]) );" << endl;
  }
}

void AlVin::print_vectors() const {
  cout << "Vectors found: " << endl;

  for (unsigned int i(0); i < vectorsCount; i++) {
    cout << "\te_" << (i + 1) << " = {";
    for (unsigned int j(0); j <= dimension; j++)
      cout << (j ? "," : "") << *vectors[i][j];
    cout << "}" << endl;
  }
}

vector<vector<unsigned int>> AlVin::get_coxeterMatrix() const {
  return coxeterMatrix;
}

unsigned int AlVin::get_dimension() const { return dimension; }

CoxIter *AlVin::get_ptrCI() const { return ptrCI; }

string AlVin::get_strCoxeterMatrix() const {
  string strCox;

  for (vector<vector<unsigned int>>::const_iterator itRow(
           coxeterMatrix.begin());
       itRow != coxeterMatrix.end(); ++itRow) {
    if (itRow != coxeterMatrix.begin())
      strCox += "\n";

    for (vector<unsigned int>::const_iterator itCol(itRow->begin());
         itCol != itRow->end(); ++itCol)
      strCox += (itCol == itRow->begin() ? "" : ",") + to_string(*itCol);
  }

  return strCox;
}

string AlVin::get_strQF(const string &strSeparator) const {
  string strQF;

  for (auto it(qf.begin()); it != qf.end(); it++) {
    if (it == qf.begin()) {
      unique_ptr<AlgebraicInteger> ptr((*it)->copy());
      ptr->opp();
      strQF = ptr->to_string();
    } else
      strQF += strSeparator + (*it)->to_string();
  }

  return strQF;
}

string AlVin::get_strAlgebraicIntegerType() const {
  return strAlgebraicIntegerType;
}

vector<std::vector<AlgebraicInteger *>> AlVin::get_vectors() const {
  return vectors;
}

unsigned int AlVin::get_vectorsCount() const { return vectorsCount; }

vector<AlgebraicInteger *> AlVin::get_qf() const { return qf; }

std::vector<AlgebraicInteger *> AlVin::get_possibleNorm2() const {
  return vf->get_aiPossibleNorm2();
}

const std::vector<AlgebraicInteger *> *AlVin::get_ptrPossibleNorm2() const {
  return vf->get_ptraiPossibleNorm2();
}

void AlVin::set_bComputeInvariantsPolyhedron(const bool &bValue) {
  bComputeInvariantsPolyhedron = bValue;
}

void AlVin::set_iCreateImage(const int &iValue) { iCreateImage = iValue; }

string AlVin::get_strFinalInformation() const { return strFinalInformation; }
