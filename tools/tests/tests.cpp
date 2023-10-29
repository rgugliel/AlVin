#include "tests.h"

// -------------------------------------------------------------------------------------
// Test AlVin
void Tests::test_AlVin() {
  vector<AlVin_test> VTlist;
  try {
    VTlist = test_AlVin_readQFList();

    cout << "Found test files: " << VTlist.size() << endl;

    for (auto vt : VTlist) {
      cout << "Run test with quadratic form: <";
      for(unsigned int i(0); i < vt.aiQF.size(); i++)
        cout << (i ? ", ": "") << *vt.aiQF[i];
      cout << "> ; expected number of vectors: " << vt.iVectorsCount << endl;

      if (vt.strField == "Rational") {
        vector<int> iQF;
        for (auto c : vt.aiQF)
          iQF.push_back(dynamic_cast<RationalInteger *>(c)->iVal);

        RationalInteger_AlVin *v(
            new RationalInteger_AlVin(iQF, "generic", false, false));

        v->Run(0, vt.iVectorsCount, false);
        delete v;
      } else if (vt.strField == "Quadratic") {
        QuadraticInteger::set_d(vt.iField);

        vector<QuadraticInteger> qiQF;
        for (auto c : vt.aiQF) {
          QuadraticInteger *qi(dynamic_cast<QuadraticInteger *>(c));

          qiQF.push_back(*qi);
        }

        QuadraticInteger_AlVin v(qiQF, "generic", false, false);
        v.Run(0, vt.iVectorsCount, false);
      } else if (vt.strField == "RC7") {
        vector<RCyclotomic7Integer> rciQF;
        for (auto c : vt.aiQF) {
          RCyclotomic7Integer *qi(dynamic_cast<RCyclotomic7Integer *>(c));

          rciQF.push_back(*qi);
        }

        RCyclotomic7Integer_AlVin v(rciQF, "generic", false, false);
        v.Run(0, vt.iVectorsCount, false);
      }
    }
  } catch (string strE) {
    cout << "Error: " << strE << endl;
  }

  cout << "Done" << endl;
  for (auto v : VTlist) {
    for (auto c : v.aiQF)
      delete c;
  }
}

vector<AlVin_test> Tests::test_AlVin_readQFList() {
  vector<AlVin_test> VTlist;

  ifstream fileIn("../quadratic_forms_tests.in");
  if (fileIn.fail())
    throw(string("Error while opening ../quadratic_forms_tests.in"));

  string strLine;
  while (getline(fileIn, strLine)) {
    PCRERegexp regexp, regexp2;
    PCREResult regexpRes, regexpRes2;

    string strField("Rational"), strFieldSupp, strVectorsCount("0");
    vector<string> strCoeffs;
    AlVin_test vt;

    // ----------------------------------------------------------------
    // field
    regexpRes.clear();
    if (regexp.preg_match_all("k[[:space:]]*=[[:space:]]*Q\\[[[:space:]]*sqrt[["
                              ":space:]]*([[:digit:]]+)",
                              strLine, regexpRes) != 0) {
      strFieldSupp = regexpRes[1][0];
      strField = "Quadratic";
    }

    if (regexp.preg_match_all("k[[:space:]]*=[[:space:]]*RC7", strLine,
                              regexpRes) != 0)
      strField = "RC7";

    // ----------------------------------------------------------------
    // vectors max
    regexpRes.clear();
    if (regexp.preg_match_all("vcount|vmax=([[:digit:]]+)", strLine, regexpRes))
      strVectorsCount = regexpRes[1][0];
    vt.iVectorsCount = stoi(strVectorsCount);

    // ----------------------------------------------------------------
    // coefficients
    regexpRes.clear();
    if (strField == "Rational") {
      if (regexp.preg_match_all("<([-[:digit:] ,]+)>", strLine, regexpRes))
        strCoeffs = explode(",", regexpRes[1][0]);

      vector<AlgebraicInteger *> aiQF;
      for (auto c : strCoeffs)
        aiQF.push_back(new RationalInteger(stoi(c)));

      vt.strField = "Rational";
      vt.aiQF = aiQF;

      VTlist.push_back(vt);
    } else if (strField == "Quadratic") {
      if (regexp.preg_match_all("<([-[:digit:]T ,]+)>", strLine, regexpRes))
        strCoeffs = explode(",", regexpRes[1][0]);

      vector<AlgebraicInteger *> aiQF;

      for (auto c : strCoeffs) {
        str_replace(c, " ", "");
        if (regexp.preg_match_all("^([-[:digit:]]+)$", c, regexpRes)) // integer
          aiQF.push_back(new QuadraticInteger(stoi(c), 0));
        else if (regexp.preg_match_all(
                     "([+-]{0,1})([[:digit:]]+)([+-]{0,1})([[:digit:]]*)", c,
                     regexpRes)) {
          string strN1(regexpRes[1][0] + regexpRes[2][0]),
              strN2(regexpRes[3][0]);

          strN2 += regexpRes[4][0] == "" ? "1" : regexpRes[4][0];

          aiQF.push_back(new QuadraticInteger(stoi(strN1), stoi(strN2)));
        } else if (regexp.preg_match_all("^([-[:digit:]]+)[*]?T$", c,
                                         regexpRes)) // non-integer part only
        {
          if (regexpRes[1][0] == "-")
            regexpRes[1][0] = "-1";

          aiQF.push_back(new QuadraticInteger(0, stoi(regexpRes[1][0])));
        } else
          cout << "Unknown type: " << c << endl;
      }

      vt.strField = "Quadratic";
      vt.aiQF = aiQF;
      vt.iField = stoi(strFieldSupp);
      VTlist.push_back(vt);
    } else if (strField == "RC7") {
      vector<AlgebraicInteger *> aiQF;
      unsigned int iRegexpCount(0);

      if (regexp.preg_match_all("<([-[:digit:]\\[\\]\\* ,]+)>", strLine,
                                regexpRes)) {
        int iBracketsNumber(0);
        char *ptr(const_cast<char *>(regexpRes[1][0].c_str()));
        while (*ptr) {
          if (*ptr == '[')
            iBracketsNumber++;
          else if (*ptr == ']')
            iBracketsNumber--;

          if (iBracketsNumber < 0)
            throw(string("The brackets in the quadratic form don't match"));

          if (iBracketsNumber && *ptr == ',')
            *ptr = ';';

          ++ptr;
        }

        if (iBracketsNumber != 0)
          throw(string("The brackets in the quadratic form don't match"));

        if ((iRegexpCount = regexp2.preg_match_all(
                 "([[:digit:]\\*\\-]*)\\[([[:digit:];\\-]+)\\]",
                 regexpRes[1][0], regexpRes2)) > 0) {
          for (unsigned i(0); i < iRegexpCount; i++) {
            str_replace(regexpRes2[1][i], "*", "");
            if (regexpRes2[1][i] == "-")
              regexpRes2[1][i] = "-1";

            int iCoefficient(regexpRes2[1][i] == "" ? 1
                                                    : stoi(regexpRes2[1][i]));

            vector<string> strCoeffs(explode(";", regexpRes2[2][i]));
            if (strCoeffs.size() != 3)
              throw(string("RCyclotomic: Bad integer: " + regexpRes2[0][i]));

            aiQF.push_back(new RCyclotomic7Integer(
                iCoefficient * (strCoeffs[0] == "" ? 0 : stoi(strCoeffs[0])),
                 iCoefficient * (strCoeffs[1] == "" ? 0 : stoi(strCoeffs[1])),
                 iCoefficient *
                     (strCoeffs[2] == "" ? 0 : stoi(strCoeffs[2]))));

            str_replace(regexpRes[1][0], regexpRes2[0][i], "");
          }
        }

        str_replace(regexpRes[1][0], " ", "");
        vector<string> strQF(explode(",", regexpRes[1][0]));

        for (auto c : strQF) {
          if (c != "")
            aiQF.push_back(new RCyclotomic7Integer(stoi(c)));
        }
        vt.strField = "RC7";
        vt.aiQF = aiQF;
        vt.iField = 0;
        VTlist.push_back(vt);
      }
    }
  }

  return VTlist;
}

// -------------------------------------------------------------------------------------
// Test numbers

void Tests::test_RationalsNumbers() {
  cout << "RationalInteger_iSQRTQuotientsup( )" << endl;
  RationalInteger_iSQRTQuotientsup();
  cout << "\tDone" << endl;

  cout << "RationalInteger_iSQRTQuotient( )" << endl;
  RationalInteger_iSQRTQuotient();
  cout << "\tDone" << endl;
}

void Tests::test_QuadraticNumbers() {
  cout << "QuadraticInteger_gcd( )" << endl;
  QuadraticInteger_gcd();
  cout << "\tDone" << endl;

  cout << "QuadraticInteger_iSQRT_quotient( )" << endl;
  QuadraticInteger_iSQRT_quotient();
  cout << "\tDone" << endl;

  cout << "QuadraticInteger_iSQRTsup_quotient( )" << endl;
  QuadraticInteger_iSQRTsup_quotient();
  cout << "\tDone" << endl;

  cout << "QuadraticInteger_bIsLessThanInt( )" << endl;
  QuadraticInteger_bIsLessThanInt();
  cout << "\tDone" << endl;

  cout << "QuadraticInteger_bIsGreaterThanInt( )" << endl;
  QuadraticInteger_bIsGreaterThanInt();
  cout << "\tDone" << endl;

  cout << "QuadraticInteger_bIsLessThanQuadraticInteger( )" << endl;
  QuadraticInteger_bIsLessThanQuadraticInteger();
  cout << "\tDone" << endl;

  cout << "QuadraticInteger_qiMultiplyBy( )" << endl;
  QuadraticInteger_qiMultiplyBy();
  cout << "\tDone" << endl;

  cout << "QuadraticInteger_primeFactors( )" << endl;
  QuadraticInteger_primeFactors();
  cout << "\tDone" << endl;
}

void Tests::test_QuadraticNumbersBig() {
  cout << "QuadraticIntegerBig_primeFactors( )" << endl;
  QuadraticIntegerBig_primeFactors();
  cout << "\tDone" << endl;

  cout << "QuadraticIntegerBig_gcd( )" << endl;
  QuadraticIntegerBig_gcd();
  cout << "\tDone" << endl;

  cout << "QuadraticIntegerBig_iSQRT_quotient( )" << endl;
  QuadraticIntegerBig_iSQRT_quotient();
  cout << "\tDone" << endl;

  cout << "QuadraticIntegerBig_iSQRTsup_quotient( )" << endl;
  QuadraticIntegerBig_iSQRTsup_quotient();
  cout << "\tDone" << endl;

  cout << "QuadraticIntegerBig_bIsLessThanInt( )" << endl;
  QuadraticIntegerBig_bIsLessThanInt();
  cout << "\tDone" << endl;

  cout << "QuadraticIntegerBig_bIsGreaterThanInt( )" << endl;
  QuadraticIntegerBig_bIsGreaterThanInt();
  cout << "\tDone" << endl;

  cout << "QuadraticIntegerBig_bIsLessThanQuadraticInteger( )" << endl;
  QuadraticIntegerBig_bIsLessThanQuadraticIntegerBig();
  cout << "\tDone" << endl;

  cout << "QuadraticIntegerBig_qiMultiplyBy( )" << endl;
  QuadraticIntegerBig_qiMultiplyBy();
  cout << "\tDone" << endl;
}

void Tests::test_Cyclotomic7Numbers() {
  RCyclotomic7Integer::initialize();

  cout << "Cyclotomic7Integer_gcd( )" << endl;
  Cyclotomic7Integer_gcd();
  cout << "\tDone" << endl;

  cout << "Cyclotomic7Integer_primeDecomposition( )" << endl;
  Cyclotomic7Integer_primeFactors();
  cout << "\tDone" << endl;

  cout << "Cyclotomic7Integer_primeFactors( )" << endl;
  Cyclotomic7Integer_primeFactors();
  cout << "\tDone" << endl;

  cout << "Cyclotomic7Integer_lessThanZero( )" << endl;
  Cyclotomic7Integer_lessThanZero();
  cout << "\tDone" << endl;

  cout << "Cyclotomic7Integer_lessThan( )" << endl;
  Cyclotomic7Integer_lessThan();
  cout << "\tDone" << endl;

  cout << "Cyclotomic7Integer_productDivision( )" << endl;
  Cyclotomic7Integer_productDivision();
  cout << "\tDone" << endl;
}

void Tests::Cyclotomic7Integer_gcd() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-11, 11);

  for (unsigned int i(0); i < 40000; i++) {
    RCyclotomic7Integer rci1(dis(gen), dis(gen), dis(gen));
    RCyclotomic7Integer rci2(dis(gen), dis(gen), dis(gen));

    RCyclotomic7Integer rci1_bck(rci1), rci2_bck(rci2), qiGCD(rci1);

    try {
      qiGCD.gcd(&rci2);

      if (!rci1.isDivisibleBy(&qiGCD) || !rci2.isDivisibleBy(&qiGCD))
        throw(string("Cyclotomic7Integer_gcd : Should be divisible by: " +
                     rci1_bck.to_string() + "," + rci2_bck.to_string()));

      if (!rci1.divideByIfDivisible(&qiGCD) ||
          !rci2.divideByIfDivisible(&qiGCD))
        throw(string("Cyclotomic7Integer_gcd : Should be divisible by 2: " +
                     rci1_bck.to_string() + "," + rci2_bck.to_string()));

      vector<RCyclotomic7Integer> rciPrimes1(rci1.rciPrimeFactors()),
          rciPrimes2(rci2.rciPrimeFactors()), rciRes;

      sort(rciPrimes1.begin(), rciPrimes1.end());
      sort(rciPrimes2.begin(), rciPrimes2.end());
      set_intersection(rciPrimes1.begin(), rciPrimes1.end(), rciPrimes2.begin(),
                       rciPrimes2.end(), back_inserter(rciRes));

      if (rciRes.size())
        throw(string("Cyclotomic7Integer_gcd : Too many factors: " +
                     rci1_bck.to_string() + "," + rci2_bck.to_string()));
    } catch (string strE) {
      cout << "Error: " << strE << endl;
      cout << "Operands: "
           << "\n\t" << rci1_bck << ", N=" << rci1_bck.iNorm() << "\n\t"
           << rci2_bck << ", N=" << rci2_bck.iNorm() << endl;
      exit(0);
    }
  }
}

void Tests::Cyclotomic7Integer_primeDecomposition() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-11, 11);

  unsigned int iMaxPrime(0), j;

  for (unsigned int i(0); i < 40000; i++) {
    RCyclotomic7Integer rci(dis(gen), dis(gen), dis(gen));
    RCyclotomic7Integer rciBck(rci);

    if (rci.isEqualTo(0))
      continue;

    try {
      auto rciPrimes(rci.rciPrimeDecomposition());

      for (auto it : rciPrimes) {
        for (j = 0; j < it.second; j++) {
          if (!rci.divideByIfDivisible(&it.first))
            throw(string("Cyclotomic7Integer_primeDecomposition: Should be "
                         "divisible by: " +
                         rciBck.to_string() + " , " + it.first.to_string()));
        }

        if (!rci.divideByIfDivisible(&it.first))
          throw(string("Cyclotomic7Integer_primeDecomposition: Should not be "
                       "divisible by: " +
                       rciBck.to_string() + " , " + it.first.to_string()));
      }

      if (abs(rci.iNorm()) != 1)
        cout << "Error (Cyclotomic7Integer_primeDecomposition): \n\trci: "
             << rciBck << endl;
    } catch (string strE) {
      if (strE == "RCyclotomic7Integer: Cyclotomic7Integer_primeDecomposition: "
                  "Max prime number")
        iMaxPrime++;

      cout << "Error: " << strE << "\n\trci: " << rci
           << ", norm=" << rci.iNorm() << endl;
    }
  }
}

void Tests::Cyclotomic7Integer_primeFactors() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-11, 11);

  unsigned int iMaxPrime(0);

  for (unsigned int i(0); i < 40000; i++) {
    RCyclotomic7Integer rci(dis(gen), dis(gen), dis(gen));
    RCyclotomic7Integer rciBck(rci);

    if (rci.isEqualTo(0))
      continue;

    try {
      auto rciPrimes(rci.rciPrimeFactors());
      for (auto it : rciPrimes) {
        if (!rci.divideByIfDivisible(&it))
          throw(string("Cyclotomic7Integer_primeFactors: Not divisible"));

        while (rci.divideByIfDivisible(&it))
          ;
      }

      if (abs(rci.iNorm()) != 1)
        cout << "Error (Cyclotomic7Integer_primeFactors): \n\trci: " << rciBck
             << endl;
    } catch (string strE) {
      if (strE == "RCyclotomic7Integer: rciPrimeFactors: Max prime number")
        iMaxPrime++;

      cout << "Error: " << strE << "\n\trci: " << rci
           << ", norm=" << rci.iNorm() << endl;
    }
  }
}

void Tests::Cyclotomic7Integer_lessThanZero() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  for (unsigned int i(0); i < 40000000; i++) {
    RCyclotomic7Integer rci(dis(gen), dis(gen), dis(gen));

    if ((rci.to_double() < 0) != rci.isLessThan(0)) {
      cout << "Error (Cyclotomic7Integer_lessThanZero): \n\trci: " << rci
           << endl;
      cout << "\tFloat approx: " << rci.to_double() << endl;
      cout << "\tTest zero: " << (int)rci.isLessThan(0) << endl;
      exit(0);
    }
  }
}

void Tests::Cyclotomic7Integer_lessThan() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  for (unsigned int i(0); i < 40000000; i++) {
    RCyclotomic7Integer rci1(dis(gen), dis(gen), dis(gen));
    RCyclotomic7Integer rci2(dis(gen), dis(gen), dis(gen));

    if ((rci1.to_double() < rci2.to_double()) != rci1.isLessThan(rci2)) {
      cout << "Error (Cyclotomic7Integer_lessThan): \n\trci1=" << rci1
           << ", rci2=" << rci2 << endl;
      exit(0);
    }
  }
}

void Tests::Cyclotomic7Integer_productDivision() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  RCyclotomic7Integer rci0(0, 0, 0);

  for (unsigned int i(0); i < 4000000; i++) {

    RCyclotomic7Integer rci1(dis(gen), dis(gen), dis(gen));
    RCyclotomic7Integer rci2(dis(gen), dis(gen), dis(gen));
    RCyclotomic7Integer rciProduct(rci1);

    if (rci1.isEqualTo(rci0) || rci2.isEqualTo(rci0))
      continue;

    double dProd(rci1.to_double() * rci2.to_double());

    rciProduct.multiplyBy(&rci2);

    if (rciProduct.to_double() - dProd > 5 * 1e-11) {
      cout << "Error (Cyclotomic7Integer_productDivision): \n\trci1: " << rci1
           << "\n\trci2: " << rci2 << "\n\tComputed product: " << rciProduct
           << "\n\tDouble value: " << rciProduct.to_double() << endl;

      cout << "différence: " << rciProduct.to_double() - dProd << endl;

      return;
    }

    if (!rciProduct.isDivisibleBy(&rci1))
      cout << "Error (Cyclotomic7Integer_productDivision) divisibility test 1: "
              "\n\trci1: "
           << rci1 << "\n\trci2: " << rci2 << endl;

    if (!rciProduct.isDivisibleBy(&rci2))
      cout << "Error (Cyclotomic7Integer_productDivision) divisibility test 2: "
              "\n\trci1: "
           << rci1 << "\n\trci2: " << rci2 << endl;

    rciProduct.divideBy(&rci2);
    if (!rciProduct.isEqualTo(rci1))
      cout << "Error (Cyclotomic7Integer_productDivision) division: \n\trci1: "
           << rci1 << "\n\trci2: " << rci2 << endl;
  }
}

void Tests::test_Numbers() {
  test_QuadraticNumbers();
  test_QuadraticNumbersBig();
  exit(0);

  test_Cyclotomic7Numbers();
  test_QuadraticNumbers();
  test_RationalsNumbers();
}

void Tests::RationalInteger_iSQRTQuotientsup() {
  std::random_device rd, rd2;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(1, 1E9);

  std::mt19937 gen2(rd2());
  std::uniform_int_distribution<> dis2(1, 1000);

  for (unsigned int i(0); i < 400000000; i++) {
    unsigned int iA(dis(gen)), iB(dis2(gen2));
    unsigned int iRes(sqrtSupQuotient<unsigned int>(iA, iB));

    if (iB * iRes * iRes < iA || iA <= iB * (iRes - 1) * (iRes - 1)) {
#pragma omp critical
      cout << "Erreur (RationalInteger_iSQRTQuotientsup): " << iA << " / " << iB
           << " ---> " << iRes << endl;
      return;
    }
  }
}

void Tests::RationalInteger_iSQRTQuotient() {
  std::random_device rd, rd2;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(1, 1E9);

  std::mt19937 gen2(rd2());
  std::uniform_int_distribution<> dis2(1, 1000);

  for (unsigned int i(0); i < 400000000; i++) {
    unsigned int iA(dis(gen)), iB(dis2(gen2));
    unsigned int iRes(sqrtQuotient<unsigned int>(iA, iB));

    if (iB * iRes * iRes > iA || iA >= iB * (iRes + 1) * (iRes + 1)) {
#pragma omp critical
      cout << "Erreur (RationalInteger_iSQRTQuotient): " << iA << " / " << iB
           << " ---> " << iRes << endl;
      return;
    }
  }
}

void Tests::QuadraticInteger_gcd() {
  try {
    std::random_device rd;

    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(-100, 100);

    std::mt19937 gend(rd());
    std::uniform_int_distribution<> disd(2, 20);

    std::mt19937 gensmall(rd());
    std::uniform_int_distribution<> dissmall(-40, 40);

    for (unsigned int i(0); i < 4000000; i++) {
      QuadraticInteger qi1(dissmall(gensmall), dissmall(gensmall));

      unsigned int dInit(disd(gend));
      unsigned int d(iRemoveSquareFactors(dInit));
      if (d == 1)
        dInit += 1;

      if (!QuadraticInteger::isDAdmissible(dInit))
        continue;

      qi1.set_d(dInit);

      QuadraticInteger qi2(dissmall(gensmall), dissmall(gensmall));
      QuadraticInteger qi1_bck(qi1), qi2_bck(qi2), qiGCD(qi1);

      qiGCD.gcd(&qi2);

      if (!qi1.isDivisibleBy(&qiGCD) || !qi2.isDivisibleBy(&qiGCD))
        throw(string("QuadraticInteger_qiGCD : Should be divisible by: " +
                     qi1_bck.to_string() + "," + qi2_bck.to_string()));

      if (!qi1.divideByIfDivisible(&qiGCD) || !qi2.divideByIfDivisible(&qiGCD))
        throw(string("QuadraticInteger_qiGCD : Should be divisible by 2: " +
                     qi1_bck.to_string() + "," + qi2_bck.to_string()));

      vector<QuadraticInteger> qiPrimes1(qi1.qiPrimeFactors()),
          qiPrimes2(qi2.qiPrimeFactors()), qiRes;

      sort(qiPrimes1.begin(), qiPrimes1.end());
      sort(qiPrimes2.begin(), qiPrimes2.end());
      set_intersection(qiPrimes1.begin(), qiPrimes1.end(), qiPrimes2.begin(),
                       qiPrimes2.end(), back_inserter(qiRes));

      if (qiRes.size())
        throw(string("QuadraticInteger_qiGCD : Too many factors: " +
                     qi1_bck.to_string() + "," + qi1_bck.to_string()));
    }
  } catch (string strE) {
    cout << "Error: " << strE << endl;
  }
}

void Tests::QuadraticInteger_iSQRTsup_quotient() {
  std::random_device rd, rd2;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-2000, 2000);

  std::mt19937 gen2(rd2());
  std::uniform_int_distribution<> dis2(-10, 10);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 80);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticInteger qi1(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticInteger::isDAdmissible(dInit))
      continue;

    QuadraticInteger::set_d(dInit);

    QuadraticInteger qi2(dis2(gen2), dis2(gen2));

    if (qi2.isEqualTo(0))
      continue;

    QuadraticInteger qiTestSign(qi1);
    qiTestSign.multiplyBy(&qi2);
    if (qiTestSign.isLessThan(0))
      continue;

    double dQuotient(qi1.to_double() / qi2.to_double());

    long int iValue;
    if (qi1.isDivisibleBy(&qi2)) {
      QuadraticInteger qiDiv(qi1);
      qiDiv.divideBy(&qi2);
      if (qiDiv.b == 0)
        iValue = sqrtSup((unsigned long int)qiDiv.a);
      else
        iValue = ceil(sqrt(dQuotient));
    } else
      iValue = ceil(sqrt(dQuotient));

    if (iValue != QuadraticInteger::iSQRTsup_quotient(qi1, qi2)) {
      cout << "QuadraticInteger_iSQRTsup_quotient: " << endl;
      cout << qi1 << " / " << qi2 << endl;
      cout << dQuotient << " --> " << ceil(sqrt(dQuotient)) << ", "
           << QuadraticInteger::iSQRTsup_quotient(qi1, qi2) << "\n"
           << endl;

      return;
    }
  }
}

void Tests::QuadraticInteger_iSQRT_quotient() {
  std::random_device rd, rd2;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-2000, 2000);

  std::mt19937 gen2(rd2());
  std::uniform_int_distribution<> dis2(-10, 10);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 80);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticInteger qi1(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticInteger::isDAdmissible(dInit))
      continue;

    QuadraticInteger::set_d(dInit);

    QuadraticInteger qi2(dis2(gen2), dis2(gen2));

    if (qi2.isEqualTo(0))
      continue;

    QuadraticInteger qiTestSign(qi1);
    qiTestSign.multiplyBy(&qi2);
    if (qiTestSign.isLessThan(0))
      continue;

    double dQuotient(qi1.to_double() / qi2.to_double());

    long int iValue;
    if (qi1.isDivisibleBy(&qi2)) {
      QuadraticInteger qiDiv(qi1);
      qiDiv.divideBy(&qi2);
      if (qiDiv.b == 0)
        iValue = integerSqrt((unsigned long int)qiDiv.a);
      else
        iValue = floor(sqrt(dQuotient));
    } else
      iValue = floor(sqrt(dQuotient));

    if (iValue != QuadraticInteger::iSQRT_quotient(qi1, qi2)) {
#pragma omp critical
      {
        cout << qi1 << " / " << qi2 << endl;
        cout << dQuotient << "--> " << floor(sqrt(dQuotient)) << ", "
             << QuadraticInteger::iSQRT_quotient(qi1, qi2) << "\n"
             << endl;
      }
      return;
    }
  }
}

void Tests::QuadraticInteger_primeFactors() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  std::mt19937 gensmall(rd());
  std::uniform_int_distribution<> dissmall(-40, 40);

  for (unsigned int i(0); i < 10000; i++) {
    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticInteger::isDAdmissible(dInit))
      continue;

    QuadraticInteger::set_d(dInit);

    QuadraticInteger qi(dissmall(gensmall), dissmall(gensmall));
    if (qi.isEqualTo(0))
      continue;

    auto qiPF(qi.qiPrimeFactors());

    for (auto iF : qiPF) {
      if (!qi.isDivisibleBy(&iF)) {
#pragma omp critical
        cout << "ERROR (QuadraticInteger_primeFactors): " << qi
             << " pas divisible par " << iF << endl;
        break;
      }

      while (qi.isDivisibleBy(&iF))
        qi.divideBy(&iF);
    }

    if (!qi.isInvertible()) {
#pragma omp critical
      {
        cout
            << "Erreur (QuadraticInteger_primeFactors): manque un facteur dans "
            << qi << endl;
        cout << "\tFinal: " << qi << endl;
        cout << "\tNorm: " << qi.iNorm() << endl;
      }
      return;
    }
  }
}

void Tests::QuadraticInteger_bIsGreaterThanInt() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  std::mt19937 genn(rd());
  std::uniform_int_distribution<> disn(-100, 100);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticInteger qi(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticInteger::isDAdmissible(dInit))
      continue;

    qi.set_d(dInit);

    int n(disd(genn));

    if (qi.isGreaterThan(n) != (qi.to_double() > n)) {
#pragma omp critical
      cout << "Error( QuadraticInteger_bIsGreaterThanInt ): " << qi << " / "
           << n << endl;
    }
  }
}

void Tests::QuadraticInteger_qiMultiplyBy() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  std::mt19937 gensmall(rd());
  std::uniform_int_distribution<> dissmall(-40, 40);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticInteger qi1(dissmall(gensmall), dissmall(gensmall));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticInteger::isDAdmissible(dInit))
      continue;

    qi1.set_d(dInit);

    QuadraticInteger qi2(dissmall(gensmall), dissmall(gensmall));

    double dProd(qi1.to_double() * qi2.to_double());

    qi1.multiplyBy(&qi2);
    if (qi1.to_double() - dProd > 5 * 1e-11) {
      cout << "Error (QuadraticInteger_qiMultiplyBy): " << qi1 << ", " << qi2
           << " / " << qi1.to_double() - dProd << endl;
      return;
    }
  }
}

void Tests::QuadraticInteger_bIsLessThanInt() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  std::mt19937 genn(rd());
  std::uniform_int_distribution<> disn(-100, 100);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticInteger qi(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticInteger::isDAdmissible(dInit))
      continue;

    qi.set_d(dInit);

    int n(disd(genn));

    if ((qi.isLessThan(n)) != (qi.to_double() < n)) {
#pragma omp critical
      cout << "Error (QuadraticInteger_bIsLessThanInt): " << qi << " / " << n
           << endl;
    }
  }
}

void Tests::QuadraticInteger_bIsLessThanQuadraticInteger() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  for (unsigned int i(0); i < 40000000; i++) {
    QuadraticInteger qi1(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticInteger::isDAdmissible(dInit))
      continue;

    qi1.set_d(dInit);

    QuadraticInteger qi2(dis(gen), dis(gen));

    if ((qi1.isLessThan(qi2)) != (qi1.to_double() < qi2.to_double())) {
#pragma omp critical
      cout << "Error (QuadraticInteger_bIsLessThanQuadraticInteger): ("
           << (qi1.isLessThan(qi2) ? "1," : "0,")
           << (qi1.to_double() - qi2.to_double()) << ") " << qi1 << " / " << qi2
           << endl;
    }
  }
}

void Tests::QuadraticIntegerBig_gcd() {
  try {
    std::random_device rd;

    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(-100, 100);

    std::mt19937 gend(rd());
    std::uniform_int_distribution<> disd(2, 20);

    std::mt19937 gensmall(rd());
    std::uniform_int_distribution<> dissmall(-40, 40);

    for (unsigned int i(0); i < 4000000; i++) {
      QuadraticIntegerBig qi1(dissmall(gensmall), dissmall(gensmall));

      unsigned int dInit(disd(gend));
      unsigned int d(iRemoveSquareFactors(dInit));
      if (d == 1)
        dInit += 1;

      if (!QuadraticIntegerBig::isDAdmissible(dInit))
        continue;

      qi1.set_d(dInit);

      QuadraticIntegerBig qi2(dissmall(gensmall), dissmall(gensmall));
      QuadraticIntegerBig qi1_bck(qi1), qi2_bck(qi2), qiGCD(qi1);

      qiGCD.gcd(&qi2);

      if (!qi1.isDivisibleBy(&qiGCD) || !qi2.isDivisibleBy(&qiGCD))
        throw(string("QuadraticIntegerBig_qiGCD : Should be divisible by: " +
                     qi1_bck.to_string() + "," + qi2_bck.to_string()));

      if (!qi1.divideByIfDivisible(&qiGCD) || !qi2.divideByIfDivisible(&qiGCD))
        throw(string("QuadraticIntegerBig_qiGCD : Should be divisible by 2: " +
                     qi1_bck.to_string() + "," + qi2_bck.to_string()));

      vector<QuadraticIntegerBig> qiPrimes1(qi1.qiPrimeFactors()),
          qiPrimes2(qi2.qiPrimeFactors()), qiRes;

      sort(qiPrimes1.begin(), qiPrimes1.end());
      sort(qiPrimes2.begin(), qiPrimes2.end());
      set_intersection(qiPrimes1.begin(), qiPrimes1.end(), qiPrimes2.begin(),
                       qiPrimes2.end(), back_inserter(qiRes));

      if (qiRes.size())
        throw(string("QuadraticIntegerBig_qiGCD : Too many factors: " +
                     qi1_bck.to_string() + "," + qi1_bck.to_string()));
    }
  } catch (string strE) {
    cout << "Error: " << strE << endl;
  }
}

void Tests::QuadraticIntegerBig_iSQRTsup_quotient() {
  std::random_device rd, rd2;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-2000, 2000);

  std::mt19937 gen2(rd2());
  std::uniform_int_distribution<> dis2(-10, 10);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 80);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticIntegerBig qi1(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticIntegerBig::isDAdmissible(dInit))
      continue;

    QuadraticIntegerBig::set_d(dInit);

    QuadraticIntegerBig qi2(dis2(gen2), dis2(gen2));

    if (qi2.isEqualTo(0))
      continue;

    QuadraticIntegerBig qiTestSign(qi1);
    qiTestSign.multiplyBy(&qi2);
    if (qiTestSign.isLessThan(0))
      continue;

    double dQuotient(qi1.to_double() / qi2.to_double());

    long int iValue;
    if (qi1.isDivisibleBy(&qi2)) {
      QuadraticIntegerBig qiDiv(qi1);
      qiDiv.divideBy(&qi2);
      if (qiDiv.b == 0)
        iValue = sqrtSup((unsigned long int)qiDiv.a.get_ui());
      else
        iValue = ceil(sqrt(dQuotient));
    } else
      iValue = ceil(sqrt(dQuotient));

    if (iValue != QuadraticIntegerBig::iSQRTsup_quotient(qi1, qi2)) {
      cout << "QuadraticIntegerBig_iSQRTsup_quotient: " << endl;
      cout << qi1 << " / " << qi2 << endl;
      cout << dQuotient << " --> " << ceil(sqrt(dQuotient)) << ", "
           << QuadraticIntegerBig::iSQRTsup_quotient(qi1, qi2) << "\n"
           << endl;

      return;
    }
  }
}

void Tests::QuadraticIntegerBig_iSQRT_quotient() {
  std::random_device rd, rd2;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-2000, 2000);

  std::mt19937 gen2(rd2());
  std::uniform_int_distribution<> dis2(-10, 10);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 80);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticIntegerBig qi1(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticIntegerBig::isDAdmissible(dInit))
      continue;

    QuadraticIntegerBig::set_d(dInit);

    QuadraticIntegerBig qi2(dis2(gen2), dis2(gen2));

    if (qi2.isEqualTo(0))
      continue;

    QuadraticIntegerBig qiTestSign(qi1);
    qiTestSign.multiplyBy(&qi2);
    if (qiTestSign.isLessThan(0))
      continue;

    double dQuotient(qi1.to_double() / qi2.to_double());

    long int iValue;
    if (qi1.isDivisibleBy(&qi2)) {
      QuadraticIntegerBig qiDiv(qi1);
      qiDiv.divideBy(&qi2);
      if (qiDiv.b == 0)
        iValue = sqrtSup((unsigned long int)qiDiv.a.get_ui());
      else
        iValue = floor(sqrt(dQuotient));
    } else
      iValue = floor(sqrt(dQuotient));

    if (iValue != QuadraticIntegerBig::iSQRT_quotient(qi1, qi2)) {
#pragma omp critical
      {
        cout << qi1 << " / " << qi2 << endl;
        cout << dQuotient << "--> " << floor(sqrt(dQuotient)) << ", "
             << QuadraticIntegerBig::iSQRT_quotient(qi1, qi2) << "\n"
             << endl;
      }
      return;
    }
  }
}

void Tests::QuadraticIntegerBig_primeFactors() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  std::mt19937 gensmall(rd());
  std::uniform_int_distribution<> dissmall(-40, 40);

  for (unsigned int i(0); i < 400000000; i++) {
    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticIntegerBig::isDAdmissible(dInit))
      continue;

    QuadraticIntegerBig::set_d(dInit);

    QuadraticIntegerBig qi(dissmall(gensmall), dissmall(gensmall));
    if (qi.isEqualTo(0))
      continue;

    auto qiPF(qi.qiPrimeFactors());

    for (auto iF : qiPF) {
      if (!qi.isDivisibleBy(&iF)) {
#pragma omp critical
        cout << "ERROR (QuadraticIntegerBig_primeFactors): " << qi
             << " pas divisible par " << iF << endl;
        break;
      }

      while (qi.isDivisibleBy(&iF))
        qi.divideBy(&iF);
    }

    if (!qi.isInvertible()) {
#pragma omp critical
      {
        cout << "Erreur (QuadraticIntegerBig_primeFactors): manque un facteur "
                "dans "
             << qi << endl;
        cout << "\tFinal: " << qi << endl;
        cout << "\tNorm: " << qi.iNorm() << endl;
      }
      return;
    }
  }
}

void Tests::QuadraticIntegerBig_bIsGreaterThanInt() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  std::mt19937 genn(rd());
  std::uniform_int_distribution<> disn(-100, 100);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticIntegerBig qi(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticIntegerBig::isDAdmissible(dInit))
      continue;

    qi.set_d(dInit);

    int n(disd(genn));

    if (qi.isGreaterThan(n) != (qi.to_double() > n)) {
#pragma omp critical
      cout << "Error( QuadraticIntegerBig_bIsGreaterThanInt ): " << qi << " / "
           << n << endl;
    }
  }
}

void Tests::QuadraticIntegerBig_qiMultiplyBy() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  std::mt19937 gensmall(rd());
  std::uniform_int_distribution<> dissmall(-40, 40);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticIntegerBig qi1(dissmall(gensmall), dissmall(gensmall));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticIntegerBig::isDAdmissible(dInit))
      continue;

    qi1.set_d(dInit);

    QuadraticIntegerBig qi2(dissmall(gensmall), dissmall(gensmall));

    double dProd(qi1.to_double() * qi2.to_double());

    qi1.multiplyBy(&qi2);
    if (qi1.to_double() - dProd > 5 * 1e-11) {
      cout << "Error (QuadraticIntegerBig_qiMultiplyBy): " << qi1 << ", " << qi2
           << " / " << qi1.to_double() - dProd << endl;
      return;
    }
  }
}

void Tests::QuadraticIntegerBig_bIsLessThanInt() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  std::mt19937 genn(rd());
  std::uniform_int_distribution<> disn(-100, 100);

  for (unsigned int i(0); i < 400000000; i++) {
    QuadraticIntegerBig qi(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticIntegerBig::isDAdmissible(dInit))
      continue;

    qi.set_d(dInit);

    int n(disd(genn));

    if ((qi.isLessThan(n)) != (qi.to_double() < n)) {
#pragma omp critical
      cout << "Error (QuadraticIntegerBig_bIsLessThanInt): " << qi << " / " << n
           << endl;
    }
  }
}

void Tests::QuadraticIntegerBig_bIsLessThanQuadraticIntegerBig() {
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-100, 100);

  std::mt19937 gend(rd());
  std::uniform_int_distribution<> disd(2, 50);

  for (unsigned int i(0); i < 40000000; i++) {
    QuadraticIntegerBig qi1(dis(gen), dis(gen));

    unsigned int dInit(disd(gend));
    unsigned int d(iRemoveSquareFactors(dInit));
    if (d == 1)
      dInit += 1;

    if (!QuadraticIntegerBig::isDAdmissible(dInit))
      continue;

    qi1.set_d(dInit);

    QuadraticIntegerBig qi2(dis(gen), dis(gen));

    if ((qi1.isLessThan(qi2)) != (qi1.to_double() < qi2.to_double())) {
#pragma omp critical
      cout << "Error (QuadraticIntegerBig_bIsLessThanQuadraticIntegerBig): ("
           << (qi1.isLessThan(qi2) ? "1," : "0,")
           << (qi1.to_double() - qi2.to_double()) << ") " << qi1 << " / " << qi2
           << endl;
    }
  }
}
