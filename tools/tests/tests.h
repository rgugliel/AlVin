#ifndef TESTS_H
#define TESTS_H

#include "../../../CoxIter/lib/regexp.h"
#include "../../../CoxIter/lib/string.h"
#include "../../quadraticinteger_alvin.h"
#include "../../quadraticinteger_big.h"
#include "../../rationalinteger_alvin.h"
#include "../../rcyclotomic7integer_alvin.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

struct AlVin_test {
  string strField;
  unsigned int iField;

  vector<AlgebraicInteger *> aiQF;
  unsigned int iVectorsCount;
};

class Tests {
  // ---------------------------------------------------------
  // AlVin
public:
  void test_AlVin();

private:
  vector<AlVin_test> test_AlVin_readQFList();

  // ---------------------------------------------------------
  // RationalInteger, QuadraticInteger
public:
  void test_Numbers();
  void test_RationalsNumbers();
  void test_QuadraticNumbers();
  void test_QuadraticNumbersBig();
  void test_Cyclotomic7Numbers();

private:
  void RationalInteger_iSQRTQuotient();
  void RationalInteger_iSQRTQuotientsup();

  void QuadraticInteger_bIsLessThanInt();
  void QuadraticInteger_bIsGreaterThanInt();
  void QuadraticInteger_bIsLessThanQuadraticInteger();
  void QuadraticInteger_qiMultiplyBy();
  void QuadraticInteger_iSQRT_quotient();
  void QuadraticInteger_iSQRTsup_quotient();
  void QuadraticInteger_primeFactors();
  void QuadraticInteger_gcd();

  void QuadraticIntegerBig_bIsLessThanInt();
  void QuadraticIntegerBig_bIsGreaterThanInt();
  void QuadraticIntegerBig_bIsLessThanQuadraticIntegerBig();
  void QuadraticIntegerBig_qiMultiplyBy();
  void QuadraticIntegerBig_iSQRT_quotient();
  void QuadraticIntegerBig_iSQRTsup_quotient();
  void QuadraticIntegerBig_primeFactors();
  void QuadraticIntegerBig_gcd();

  void Cyclotomic7Integer_productDivision();
  void Cyclotomic7Integer_lessThanZero();
  void Cyclotomic7Integer_lessThan();
  void Cyclotomic7Integer_primeFactors();
  void Cyclotomic7Integer_primeDecomposition();
  void Cyclotomic7Integer_gcd();
};

#endif // TESTS_H
