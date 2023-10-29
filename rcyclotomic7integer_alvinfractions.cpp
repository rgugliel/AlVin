#include "rcyclotomic7integer_alvinfractions.h"

using namespace RCyclotomic7Integer_constants;

RCyclotomic7Integer_VFs::RCyclotomic7Integer_VFs(
    vector<AlgebraicInteger *> aiPossibleNorms2,
    const RCyclotomic7Integer &rciAlpha0)
    : rciAlpha0(rciAlpha0), rciMinusAlpha0(rciAlpha0),
      AlVinFractions(aiPossibleNorms2),
      rciPossibleNorms2_max(dynamic_cast<RCyclotomic7Integer *>(
          aiPossibleNorms2[aiPossibleNorms2.size() - 1])) {
  rciMinusAlpha0.multiplyBy(-1);
}

RCyclotomic7Integer_VFs::~RCyclotomic7Integer_VFs() {}

void RCyclotomic7Integer_VFs::computeNextAlVinFractions() {
  // ----------------------------------------
  // Frequent values
  static interval gaol_minus1(interval(-1));
  static interval gaol_2(interval(2));
  static interval gaol_3(interval(3));
  static interval gaol_7(interval(7));

  // ----------------------------------------
  // Bounds
  unsigned int iMin(iLastMaximum),
      iSizeTemp((iMin + iBatchSize) * (iMin + iBatchSize) - iLastMaximum);
  iLastMaximum += min(iSizeTemp, (unsigned int)10000);

  long int iTemp, iSqrtTemp;
  interval gaol_b, gaol_a;

  // We generate fractions f such that: iMin < f <= iLastMaximum
  for (auto aiE : aiPossibleNorms2) {
    // We write the numerator: a * l1 + b * l2 + c * l3
    RCyclotomic7Integer *rciE(dynamic_cast<RCyclotomic7Integer *>(aiE));

    interval gaol_Mprime(sqrt(rciE->to_interval() * interval(iLastMaximum)));
    interval gaol_mprime(sqrt(rciE->to_interval() * interval(iMin)));

    bool bIsMprimeLongInt(false), bIsmprimeLongInt(false);
    long int iMprime(0),
        imprime(0); // Used only if bIsMprimeLongInt/bIsmprimeLongInt is true

    if (rciE->isLongInt()) {
      // If m' is an integer
      iTemp = -rciE->iC[0].get_si() * iMin;

      iSqrtTemp =
          integerSqrt<unsigned long int>(iTemp); // rciE & iLastMaximum are positive
      if (iSqrtTemp * iSqrtTemp == iTemp)  // if iTemp is a square
      {
        bIsmprimeLongInt = true;
        imprime = iSqrtTemp;
      }

      // If M' is an integer
      iTemp = -rciE->iC[0].get_si() * iLastMaximum;

      iSqrtTemp =
          integerSqrt<unsigned long int>(iTemp); // rciE & iLastMaximum are positive
      if (iSqrtTemp * iSqrtTemp == iTemp)  // if iTemp is a square
      {
        bIsMprimeLongInt = true;
        iMprime = iSqrtTemp;
      }
    }

    // ------------------------------------
    // Bounds for b
    // First term: (2 * l1 + 2 * l2 + 3 * l3 ) * R2
    rciE->conjugate(2);
    rciMinusAlpha0.conjugate(2);
    interval gaol_R2(sqrt(rciE->to_interval() / rciMinusAlpha0.to_interval()));
    interval gaol_temp((gaol_2 * (gaol_l1 + gaol_l2) + gaol_3 * gaol_l3) *
                       gaol_R2);

    // Second term: (l1 + 4 * l2 + 2 * l3 ) ) * R3
    rciE->conjugate(2);
    rciMinusAlpha0.conjugate(2);
    interval gaol_R3(sqrt(rciE->to_interval() / rciMinusAlpha0.to_interval()));
    gaol_temp += (gaol_l1 + gaol_2 * (gaol_2 * gaol_l2 + gaol_l3)) * gaol_R3;

    // Third term (min):  M' * ( 2 * l1 + 3 * l2 + 2 * l3 )
    interval gaol_min(gaol_temp + gaol_Mprime * (gaol_2 * (gaol_l1 + gaol_l3) +
                                                 gaol_3 * gaol_l2));

    // Third term (max):
    interval gaol_max(gaol_temp * gaol_minus1 +
                      gaol_mprime *
                          (gaol_2 * (gaol_l1 + gaol_l3) + gaol_3 * gaol_l2));

    rciMinusAlpha0.conjugate(2); // Get back to the original value

    gaol_min = ceil(gaol_min / gaol_7);
    gaol_max = floor(gaol_max / gaol_7);

    if (!gaol_min.is_an_int() || !gaol_max.is_an_int())
      throw(string("RCyclotomic7Integer_VFs::computeNextAlVinFractions() : Not "
                   "an int (bounds for b)"));

    long int bMax(floor(gaol_max.left()));

    for (long int b(floor(gaol_min.left())); b <= bMax; b++) {
      gaol_b = interval(b);

      // ------------------------------------
      // Bounds for a
      gaol_min =
          (gaol_minus1 * (gaol_2 * (gaol_l1 + gaol_2 * gaol_l3) + gaol_l2));
      gaol_max = gaol_min;

      gaol_min *= ((gaol_2 * (gaol_l1 + gaol_l3) + gaol_l2 * gaol_3) * gaol_b +
                   gaol_l3 * gaol_R2 - gaol_R3 - gaol_Mprime);
      gaol_max *= ((gaol_2 * (gaol_l1 + gaol_l3) + gaol_l2 * gaol_3) * gaol_b -
                   gaol_l3 * gaol_R2 + gaol_R3 - gaol_mprime);

      gaol_min = ceil(gaol_min / gaol_7);
      gaol_max = floor(gaol_max / gaol_7);

      if (!gaol_min.is_an_int() || !gaol_max.is_an_int())
        throw(string("RCyclotomic7Integer_VFs::computeNextAlVinFractions() : "
                     "Not an int (bounds for a)"));

      long int aMax(floor(gaol_max.left()));

      for (long int a(floor(gaol_min.left())); a <= aMax; a++) {
        gaol_a = interval(a);

        long cMin, cMax;

        // ------------------------------------
        // Bounds for c
        // Inf
        if (bIsMprimeLongInt && a == b && b == -iMprime)
          cMin = a;
        else {
          gaol_min =
              (gaol_Mprime - gaol_a * gaol_l1 - gaol_b * gaol_l2) / gaol_l3;
          gaol_min = ceil(gaol_min);

          if (!gaol_min.is_an_int())
            throw(string("RCyclotomic7Integer_VFs::computeNextAlVinFractions() "
                         ": Not an int (bound inf for c)"));

          cMin = floor(gaol_min.left());
        }

        // Sup
        if (bIsmprimeLongInt && a == b && b == -imprime)
          cMax = a - 1;
        else {
          gaol_max =
              (gaol_mprime - gaol_a * gaol_l1 - gaol_b * gaol_l2) / gaol_l3;
          gaol_max = floor(gaol_max);

          if (!gaol_max.is_an_int())
            throw(string("RCyclotomic7Integer_VFs::computeNextAlVinFractions() "
                         ": Not an int (bound sup for c)"));

          cMax = floor(gaol_max.left());
        }

        for (long int c(cMin); c <= cMax; c++) {
          RCyclotomic7Integer rci(a, b, c);
          interval gaol_rci(rci.to_interval());

          if (gaol_rci.certainly_leq(gaol_mprime))
            continue;
          if (!gaol_mprime.certainly_le(
                  gaol_rci)) // If we cannot decide if m' < k0
            throw(string("RCyclotomic7Integer_VFs::computeNextAlVinFractions():"
                         " Impossible to compare"));

          rci.conjugate(2);
          gaol_temp = abs(rci.to_interval());
          if (gaol_temp.certainly_ge(gaol_R2))
            continue;
          if (!gaol_temp.certainly_leq(
                  gaol_R2)) // If we cannot decide if | sigma_2( k0 ) | < S_2
            throw(string("RCyclotomic7Integer_VFs::computeNextAlVinFractions():"
                         " Impossible to compare"));

          rci.conjugate(2);
          gaol_temp = abs(rci.to_interval());
          if (gaol_temp.certainly_ge(gaol_R3))
            continue;
          if (!gaol_temp.certainly_leq(
                  gaol_R3)) // If we cannot decide if | sigma_3( k0 ) | < S_3
            throw(string("RCyclotomic7Integer_VFs::computeNextAlVinFractions():"
                         " Impossible to compare"));

          rci.conjugate(2); // We restore the correct value

          // ------------------------------------------------
          // Add the fraction
          RCyclotomic7Integer *rciNum(new RCyclotomic7Integer(rci));

          rciNum->multiplyBy(&rci);
          rciNum->multiplyBy(rciPossibleNorms2_max);
          rciNum->divideBy(rciE);

          AlVinFraction *vf =
              new AlVinFraction(new RCyclotomic7Integer(rci),
                                new RCyclotomic7Integer(*rciE), rciNum);
          alvinfractions.insert(lower_bound(alvinfractions.begin(),
                                            alvinfractions.end(), vf,
                                            isLessThanPtrAlVinFraction),
                                vf);
        }
      }
    }
  }

  alvinfractions_it = alvinfractions.begin();
}
