#include "app.h"
#include "CoxIter/growthrate.h"

App::App()
    : bCheckNR(false), bCheckNREquations(false), bComputeInvariantsQF(false),
      bComputeInvariantsPolyhedron(false), bDebug(false), bPrintHelp(false),
      iCreateImage(-1), iMaxVectors(0), iMinVectors(0), iFieldSupp(0),
      strField("rationals"), strOuputMathematicalFormat("generic") {}

App::~App() {
  for (auto c : aiQF)
    delete c;
}

void App::readMainParameters(int argc, char **argv) {
  string strTemp;
  PCRERegexp regexp, regexp2;
  PCREResult regexpRes, regexpRes2;
#ifdef _RC7AVAILABLE_
  int iRegexpCount;
#endif

  // --------------------------------------------------
  // Complete parameters
  string strParams;
  for (int i(0); i < argc; ++i)
    strParams += " " + std::string(argv[i]);

  // --------------------------------------------------
  // Debug?
  if (strParams.find("-h") != string::npos ||
      strParams.find("-help") != string::npos || argc == 1) {
    bPrintHelp = true;
    return;
  }

  // --------------------------------------------------
  // Field of definition
  if (regexp.preg_match_all(
          "([-]*[[:space:]]*k|-[[:space:]]*field)[[:space:]]*[=]*[[[:space:]Q]*"
          "[\\[\\(]?[[:space:]]*sqrt[[:space:]\\[\\(]*([0-9]{1,5})",
          strParams, regexpRes, PCRE_CASELESS)) {
    iFieldSupp = stoi(regexpRes[2][0]);

    if (!QuadraticInteger::bIsDAdmissible(iFieldSupp))
      throw(string("Field of definition: quadratic field: the following value "
                   "is not admissible: " +
                   regexpRes[2][0]));

    strField = "quadratic";
  }
#ifdef _RC7AVAILABLE_
  else if (regexp.preg_match_all("(-k|-[[:space:]]*k|-[[:space:]]*field)[[:"
                                 "space:]]*[=]*[[:space:]]*RC7",
                                 strParams, regexpRes, PCRE_CASELESS)) {
    strField = "rc7";
  }
#endif

  // --------------------------------------------------
  // Debug?
  if (strParams.find("-debug") != string::npos) {
    str_replace(strParams, "-debug", "");
    bDebug = true;
  }

  // --------------------------------------------------
  // Image
  if (strParams.find("-noimage") != string::npos) {
    str_replace(strParams, "-noimage", "");
    iCreateImage = 0;
  }
  if (strParams.find("-image") != string::npos) {
    str_replace(strParams, "-image", "");
    iCreateImage = 1;
  }

  // --------------------------------------------------
  // Minimal number of vectors to compute
  if (regexp.preg_match_all("-minv[[:space:]]+([[:digit:]]{1,4})", strParams,
                            regexpRes)) {
    str_replace(strParams, regexpRes[0][0], "");
    iMinVectors = stoi(regexpRes[1][0]);
  }

  // --------------------------------------------------
  // Maximal number of vectors to compute
  if (regexp.preg_match_all("-maxv[[:space:]]+([[:digit:]]{1,4})", strParams,
                            regexpRes)) {
    str_replace(strParams, regexpRes[0][0], "");
    iMaxVectors = stoi(regexpRes[1][0]);
  }

  // --------------------------------------------------
  // Non reflective with symmetries
  if (regexp.preg_match_all(
          "-nr[[:space:]]?\\[[[:space:]]?([[:digit:]]+)[[:space:]]?,[[:space:]]"
          "?([[:digit:]]+)\\][[:space:]]?",
          strParams, regexpRes, PCRE_CASELESS)) {
    bCheckNR = true;
    iNRMin = min(stoi(regexpRes[1][0]), stoi(regexpRes[2][0]));
    iNRMax = max(stoi(regexpRes[1][0]), stoi(regexpRes[2][0]));
    str_replace(strParams, regexpRes[0][0], "");
  }

  // --------------------------------------------------
  // Non reflective with equations
  if (strParams.find("-nrequations") != string::npos) {
    str_replace(strParams, "-nrequations", "");
    bCheckNREquations = true;
  }

  // --------------------------------------------------
  // Output format
  if (regexp.preg_match_all(
          "-oformat[[:space:]=]?(mathematica|generic|latex|pari)", strParams,
          regexpRes, PCRE_CASELESS)) {
    std::transform(regexpRes[1][0].begin(), regexpRes[1][0].end(),
                   regexpRes[1][0].begin(), ::tolower);

    strOuputMathematicalFormat = regexpRes[1][0];

    str_replace(strParams, regexpRes[0][0], "");
  }

  // --------------------------------------------------
  // Compute invariants of the polyhedron
  if (strParams.find("-ip") != string::npos ||
      strParams.find("-ipolyhedron") != string::npos ||
      strParams.find("-invariantpolyhedron") != string::npos ||
      strParams.find("-invariantspolyhedron") != string::npos) {
    str_replace(strParams, "-ip", "");
    str_replace(strParams, "-ipolyhedron", "");
    str_replace(strParams, "-invariantpolyhedron", "");
    str_replace(strParams, "-invariantspolyhedron", "");

    bComputeInvariantsPolyhedron = true;
  }

  // --------------------------------------------------
  // Invariants of quadratic form
  if (strParams.find("-iqf") != string::npos && strField == "rationals") {
    str_replace(strParams, "-iqf", "");
    bComputeInvariantsQF = true;
  }

  regexpRes.clear();

  // --------------------------------------------------
  // Quadratic form
  if (strField == "rationals") {
    if (regexp.preg_match_all("(-f|-qf|- f|- qf)([[:digit:]-, ]+)", strParams,
                              regexpRes)) {
      strTemp = regexpRes[2][0];
      str_replace(strTemp, " ", "");

      vector<string> strQF(explode(",", strTemp));

      try {
        for (auto strCo : strQF)
          aiQF.push_back(new RationalInteger(stoi(strCo)));
      } catch (std::exception &ex) {
        cout << "Quadratic form: unknown coefficient: " << strTemp << endl;
        exit(0);
      }
    }
  } else if (strField == "quadratic") {
    if (regexp.preg_match_all(
            "(-[[:space:]]*f|-[[:space:]]*qf)([[:digit:],Tt\\+\\-\\* ]+)",
            strParams, regexpRes)) {
      strTemp = regexpRes[2][0];
      str_replace(strTemp, " ", "");

      vector<string> strQF(explode(",", strTemp));

      for (auto c : strQF) {
        if (regexp.preg_match_all("^([-[:digit:]]+)$", c, regexpRes)) // integer
          aiQF.push_back(new QuadraticInteger(stoi(c), 0));
        else if (regexp.preg_match_all("([+-]{0,1})([[:digit:]]+)([+-]{1,1})([["
                                       ":digit:]]*)[\\*]*T",
                                       c, regexpRes)) {
          string strN1(regexpRes[1][0] + regexpRes[2][0]),
              strN2(regexpRes[3][0]);

          strN2 += regexpRes[4][0] == "" ? "1" : regexpRes[4][0];

          aiQF.push_back(new QuadraticInteger(stoi(strN1), stoi(strN2)));
        } else if (regexp.preg_match_all("([+-]{0,1})([[:digit:]]*)[\\*]*T([+-]"
                                         "{1,1})([[:digit:]]*)",
                                         c, regexpRes)) {
          string strN1(regexpRes[1][0]), strN2(regexpRes[3][0]);

          strN1 += regexpRes[2][0] == "" ? "1" : regexpRes[2][0];
          strN2 += regexpRes[4][0] == "" ? "1" : regexpRes[4][0];

          aiQF.push_back(new QuadraticInteger(stoi(strN2), stoi(strN1)));
        } else if (regexp.preg_match_all("^([-[:digit:]]+)[*]?(T|t)$", c,
                                         regexpRes)) // non-integer part only
        {
          if (regexpRes[1][0] == "-")
            regexpRes[1][0] = "-1";

          aiQF.push_back(new QuadraticInteger(0, stoi(regexpRes[1][0])));
        } else
          throw(string("Quadratic form: unknown coefficient: " + c));
      }
    }
  }
#ifdef _RC7AVAILABLE_
  else if (strField == "rc7") {
    if (regexp.preg_match_all(
            "(-f|-qf|- f|- qf)([[:digit:],\\+\\-\\* \\[\\]\\(\\)]+)", strParams,
            regexpRes)) {
      str_replace(regexpRes[0][0], " ", "");
      str_replace(regexpRes[0][0], "-f", "");
      str_replace(regexpRes[0][0], "-qf", "");

      regexpRes2.clear();

      // -----------------------------------------------
      // We replace , inside [] by ;
      int iBracketsNumber(0);
      char *ptr(const_cast<char *>(regexpRes[0][0].c_str()));
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
               "([[:digit:]\\*\\-]*)\\[([[:digit:];\\-]+)\\]", regexpRes[0][0],
               regexpRes2)) > 0) {
        for (int i(0); i < iRegexpCount; i++) {
          str_replace(regexpRes2[1][i], "*", "");
          if (regexpRes2[1][i] == "-")
            regexpRes2[1][i] = "-1";

          int iCoefficient(regexpRes2[1][i] == "" ? 1 : stoi(regexpRes2[1][i]));

          vector<string> strCoeffs(explode(";", regexpRes2[2][i]));
          if (strCoeffs.size() != 3)
            throw(string("RCyclotomic: Bad integer: " + regexpRes2[0][i]));

          aiQF.push_back(new RCyclotomic7Integer(
              {iCoefficient * (strCoeffs[0] == "" ? 0 : stoi(strCoeffs[0])),
               iCoefficient * (strCoeffs[1] == "" ? 0 : stoi(strCoeffs[1])),
               iCoefficient * (strCoeffs[2] == "" ? 0 : stoi(strCoeffs[2]))}));

          str_replace(regexpRes[0][0], regexpRes2[0][i], "");
        }
      }

      vector<string> strQF(explode(",", regexpRes[0][0]));

      for (auto c : strQF) {
        if (c != "")
          aiQF.push_back(new RCyclotomic7Integer(stoi(c)));
      }
    }
  }
#endif

#ifndef _DOT_PROGRAM_FOUND_
  iCreateImage = 0;
#endif
}

AlVin *App::instanciateAlVin() {
  AlVin *v(nullptr);

  if (strField == "rationals") {
    vector<int> iQF;
    for (auto c : aiQF)
      iQF.push_back(dynamic_cast<RationalInteger *>(c)->iVal);

    v = new RationalInteger_AlVin(iQF, strOuputMathematicalFormat, true,
                                  bDebug);
  } else if (strField == "quadratic") {
    QuadraticInteger::set_d(iFieldSupp);

    vector<QuadraticInteger> qiQF;
    for (auto c : aiQF) {
      QuadraticInteger *qi(dynamic_cast<QuadraticInteger *>(c));

      qiQF.push_back(*qi);
    }

    v = new QuadraticInteger_AlVin(qiQF, strOuputMathematicalFormat, true,
                                   bDebug);
  }
#ifdef _RC7AVAILABLE_
  else if (strField == "rc7") {
    vector<RCyclotomic7Integer> rciQF;
    for (auto c : aiQF) {
      RCyclotomic7Integer *rci(dynamic_cast<RCyclotomic7Integer *>(c));

      rciQF.push_back(*rci);
    }

    v = new RCyclotomic7Integer_AlVin(rciQF, strOuputMathematicalFormat, true,
                                      bDebug);
  }
#endif

  v->set_bComputeInvariantsPolyhedron(bComputeInvariantsPolyhedron);
  v->set_iCreateImage(iCreateImage);

  return v;
}

NotReflective *App::instanciateNotReflectiveEquations(AlVin *v) {
  NotReflective *nr(nullptr);

  if (strField == "rationals") {
    nr = new RationalInteger_NotReflective(v);
  } else
    throw(string("Not implemented for this field"));

  return nr;
}

InfiniteNSymetries *App::instanciateInfiniteNSymetries(AlVin *v) {
  InfiniteNSymetries *ins(nullptr);

  if (!iMaxVectors)
    throw(string("Option -nr mut be used with -maxv"));

  if (iMaxVectors < v->get_iDimension())
    throw(string(
        "The number of vectors should be at least equal to the dimension"));

  iMaxVectors = max(iNRMax, v->get_iDimension() + 1);
  iNRMin = max(iNRMin, v->get_iDimension() + 1);
  iNRMax = min(iMaxVectors, iNRMax);

  if (strField == "rationals") {
    ins = new RationalInteger_InfiniteNSymetries(v);
  } else if (strField == "quadratic") {
    ins = new QuadraticInteger_InfiniteNSymetries(v);
  }
#ifdef _RC7AVAILABLE_
  else if (strField == "rc7") {
    ins = new RCyclotomic7Integer_InfiniteNSymetries(v);
  }
#endif
  else
    throw(string("Not implemented for this field"));

  return ins;
}

void App::Run() {
  if (bPrintHelp) {
    printHelp();
    return;
  }

  if (!aiQF.size())
    throw(string("No quadratic form given"));

  try {
    AlVin *v(instanciateAlVin());

    NotReflective *nrEquations(nullptr);
    InfiniteNSymetries *ins(nullptr);

    auto r01 = RCyclotomic7Integer({1, 1, 1});
    auto r12 = RCyclotomic7Integer({3, 1, 2});
    auto r20 = RCyclotomic7Integer({2, 0, 1});
    auto r21 = RCyclotomic7Integer({2, 0, 1});

    auto a21 = RCyclotomic7Integer({0, 0, -1});

    auto a30 = RCyclotomic7Integer({0, -1, -2});
    auto a31 = RCyclotomic7Integer({0, 0, -1});
    auto a32 = RCyclotomic7Integer({1, 0, 0});
    const vector<AlgebraicInteger*> a3 = {&a30, &a31, &a32};

    auto r30 = RCyclotomic7Integer({-3, -1, -3});
    auto r31 = RCyclotomic7Integer(0);
    auto r32 = RCyclotomic7Integer({-4, -1, -3});
    const vector<AlgebraicInteger*> r3 = {&r30, &r31, &r32};

    auto qCoeff = RCyclotomic7Integer({-2, -5, -7});
    auto qRNormStrange = RCyclotomic7Integer({-3, -2, -2});

    auto display = [](AlVin *v, const vector<AlgebraicInteger*> vec){
      auto prod = dynamic_cast<RCyclotomic7Integer*>(v->aiBilinearProduct(vec, vec));
      auto factors = prod->rciPrimeDecomposition();

      cout << "\n";
      cout << "Norm^2: " << *prod << endl;
      cout << "Invertible: " << prod->isInvertible() << endl;
      cout << "Prime factors: " << factors.size() << endl;

      for(const auto& f: factors)
        cout << "\t" << f.first << " -> " << f.second << endl;

      cout << *vec[0] << "/" << *prod << endl;

      delete prod;
    };

    cout << "Coeff #prime: " << qCoeff.rciPrimeDecomposition().size() << endl;
    cout << "Related: " << qCoeff.isDivisibleBy(&qRNormStrange) << endl;

    display(v, a3);
    display(v, r3);

    chrono::time_point<std::chrono::system_clock> timeStart(
        chrono::system_clock::now());
    if (!v->Run(iMinVectors, iMaxVectors)) {
      cout << "\nThe algorithm did not terminate; the polyhedron may be of "
              "infinite volume\n"
           << endl;

      if (bCheckNREquations) {
        if (iMaxVectors) {
          nrEquations = instanciateNotReflectiveEquations(v);
          nrEquations->Run();
          delete nrEquations;
        } else
          cout << "Option -nrequations mut be used with -maxv" << endl;
      }

      if (bCheckNR) {
        try {
          ins = instanciateInfiniteNSymetries(v);
          cout << "Checking if the form is non-reflective..." << endl;
          if (ins->Run(iNRMin, iNRMax)) {
            cout << "\tThe form is non-reflective" << endl;

            if (bDebug) {
              vector<GraphInvolution> grui(ins->get_usefulInvolutions());

              cout << "\tList of used involutions:" << endl;

              for (auto inv : grui) {
                unsigned int iSize(inv.iPermutation.size());

                if (strOuputMathematicalFormat == "latex") {
                  for (unsigned int i(0); i < iSize; i++) {
                    if (inv.iPermutation[i] == i)
                      cout << (i ? ", " : "\t\t") << "e_"
                           << (inv.iVertices[i] < 9 ? "" : "{")
                           << (inv.iVertices[i] + 1)
                           << (inv.iVertices[i] < 9 ? "" : "}");
                    else if (i < inv.iPermutation[i])
                      cout << (i ? ", " : "\t\t") << "e_"
                           << (inv.iVertices[i] < 9 ? "" : "{")
                           << (inv.iVertices[i] + 1)
                           << (inv.iVertices[i] < 9 ? "" : "}")
                           << " \\leftrightarrow e_"
                           << (inv.iVertices[inv.iPermutation[i]] < 9 ? ""
                                                                      : "{")
                           << (inv.iVertices[inv.iPermutation[i]] + 1)
                           << (inv.iVertices[inv.iPermutation[i]] < 9 ? ""
                                                                      : "}");
                  }
                } else {
                  for (unsigned int i(0); i < iSize; i++) {
                    if (inv.iPermutation[i] == i)
                      cout << (i ? ", " : "\t\t") << "e"
                           << (inv.iVertices[i] + 1);
                    else if (i < inv.iPermutation[i])
                      cout << (i ? ", " : "\t\t") << "e"
                           << (inv.iVertices[i] + 1) << " <-> e"
                           << (inv.iVertices[inv.iPermutation[i]] + 1);
                  }
                }
                cout << endl;
              }
            }
          } else {
            cout << "\tCannot decide" << endl;

            if (bDebug) {
              cout << "\tFor the found involutions, a basis of the fixed "
                      "points space is:"
                   << endl;
              ins->print_basisFixedPoints("\t\t");
            }
          }

          delete ins;
        } catch (string strE) {
          cout << strE << endl;
        }
      }
    } else {
      if (bComputeInvariantsPolyhedron) {
        CoxIter *ci(v->get_ptrCI());
        ci->set_ouputMathematicalFormat(strOuputMathematicalFormat);
        ci->computeEulerCharacteristicFVector();

        unsigned int iDimension(v->get_iDimension());

        cout << "\n---------------------------------\nInformation about the "
                "polyhedron:\n---------------------------------"
             << endl;
        cout << "Euler characteristic: " << ci->get_brEulerCaracteristic()
             << endl;

        // ----------------------------------------------
        // Covolume
        if (iDimension % 2 == 0) {
          cout << "Volume: ";
          MPZ_rational cov((iDimension / 2) % 2 ? -1 : 1);
          for (unsigned int i(1); i <= iDimension; i++) {
            cov *= 2;
            cov /= i;
            if (i <= (iDimension / 2))
              cov *= i;
          }
          cout << "pi^" << (iDimension / 2) << " * "
               << cov * ci->get_brEulerCaracteristic() << endl;
        }

        // ----------------------------------------------
        // f-vector
        vector<unsigned int> iFVector(ci->get_fVector());
        cout << "f-vector: (";
        for (unsigned int i(0); i <= iDimension; i++)
          cout << (i ? ", " : "") << iFVector[i];
        cout << ")" << endl;

        cout << "Number of vertices at infinity: "
             << ci->get_verticesAtInfinityCount() << endl;

        // ---------------------------------------------
        // Number of faces with x vertices
        /*
        cout << "Number of facets with x vertices (x>=3): ";
        vector< unsigned int > iFacesAppears( ci->get_iVerticesCount(), 0 );
        auto pG( ci->get_ptr_graphsProducts() );
        for( auto it : (*pG)[2] )
        {
                for( auto v : it.get_iVertices() )
                        iFacesAppears[v]++;
        }

        for( unsigned int i(3); i < iFVector[0]; i++ )
        {
                unsigned int iTemp(0);
                for( auto it : iFacesAppears )
                {
                        if( it == i )
                                iTemp++;
                }

                cout << iTemp << ", ";
        }
        cout << endl;*/

        // ---------------------------------------------
        // Growth series and growth rate
        cout << "\nGrowth series: " << endl;
        ci->printGrowthSeries();
        cout << endl;

        GrowthRate gr;
        GrowthRate_Result grr(
            gr.grrComputations(ci->get_growthSeries_denominator()));
        if (grr.isComputed && ci->get_isGrowthSeriesReduced()) {
          cout << "\nGrowth rate: " << grr.growthRate << endl;
          cout << "\tPerron number: "
               << (grr.perron < 0 ? "?" : (grr.perron > 0 ? "yes" : "no"))
               << endl;
          cout << "\tPisot number: "
               << (grr.pisot < 0 ? "?" : (grr.pisot > 0 ? "yes" : "no"))
               << endl;
          cout << "\tSalem number: "
               << (grr.salem < 0 ? "?" : (grr.salem > 0 ? "yes" : "no"))
               << endl;
        }
      }
    }

    // Here we have K=Q
    if (bComputeInvariantsQF) {
      cout << "\nCommensurability invariant:" << endl;
      vector<int> iQF;
      for (auto ai : aiQF)
        iQF.push_back(dynamic_cast<RationalInteger *>(ai)->get_iValue());

      InvariantsQF invqf(iQF);
      cout << "\t" << invqf.get_strInvariant() << endl;
    }

    string strFinalInformation(v->get_strFinalInformation());
    if (strFinalInformation != "")
      cout << "\n" << strFinalInformation << endl;

    cout << "\nComputation time: "
         << chrono::duration<double, milli>(chrono::system_clock::now() -
                                            timeStart)
                    .count() /
                1000
         << "s" << endl;

    if (v != nullptr)
      delete v;
  } catch (string &strE) {
    cout << "Error: " << strE << endl;
    return;
  }
}

void App::printHelp() const {
  cout << "           _ __      __ _        \n"
          "    /\\    | |\\ \\    / /(_)       \n"
          "   /  \\   | | \\ \\  / /  _  _ ___  \n"
          "  / /\\ \\  | |  \\ \\/ /  | || '_  \\ \n"
          " / ____ \\ | |   \\  /   | || | | |\n"
          "/_/    \\_\\|_|    \\/    |_||_| |_|\n"
       << endl;

  cout << "AlVin is an implementation of the Vinberg algorithm\n"
          "for the fields Q, Q[sqrt(d)] and Q(cos(2*pi/7))\n\n"
          "One basic example is the following: \n"
          "\t./alvin -qf -1,2,2,3\n"
          "\tHere, we applied the algorithm to the diagonal quadratic"
          "\n\tform <-1,2,2,3>\n\n"
          "If we want to work over the quadratic field Q[sqrt(2)], we can use\n"
          "\t./alvin -k=Q[sqrt 2] -qf -1-T,1,1,1,1\n"
          "\twhere T denotes the generator of the ring of integers\n\n"
          "Another example over Q[sqrt(3)] is given by\n"
          "\t./alvin -k=Q[sqrt 3] -qf -3-2T,1,1,1,1 -ip\n"
          "\tthe -ip parameter asks AlVin to compute the Invariants of the\n"
          "\tcorresponding Polyhedron\n\n"
          "The complete documentation is available here:\n"
          "\thttps://rgugliel.github.io/AlVin/\n"
       << endl;
}
