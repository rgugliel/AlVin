#include <iomanip>
#include <iostream>

#include "tests.h"

#include <gmpxx.h>

using namespace std;

int main(int argc, char **argv) {
  /*
  QuadraticIntegerBig::set_d(5);
  QuadraticInteger::set_d(5);

  QuadraticInteger qis1( 720, 1440 ), qis2( 2, 4 );
  QuadraticIntegerBig qi1( 720, 1440 ), qi2( 2, 4 );

  cout << qis1 << ", " << qis2 << endl;
  QuadraticInteger::iSQRTsup_quotient( qis1, qis2 );

  cout << "---------------------\n" << endl;

  cout << qi1 << ", " << qi2 << endl;
  QuadraticIntegerBig::iSQRTsup_quotient( qi1, qi2 );


  return 0;*/

  /*
  RCyclotomic7Integer rci( {-61, 57, -97} );
  cout << rci.to_double() << endl;
  return 0;*/

  Tests t;
  /*
  cout << "-----------------------------------------------" << endl;
  cout << "Test numbers" << endl;
  cout << "-----------------------------------------------" << endl;
  t.test_Numbers( );*/

  t.test_AlVin();
  return 0;

  /*
  cout << "-----------------------------------------------" << endl;
  cout << "Test numbers" << endl;
  cout << "-----------------------------------------------" << endl;
  t.test_Numbers( );*/

  return 0;
}
