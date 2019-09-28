#include <chrono>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#include "app.h"

int main(int argc, char **argv) {
  try {
    App app;

    // Reading the parameters
    app.readMainParameters(argc, argv);

    app.Run();

    return 0;
  } catch (string &strE) {
    cout << "\n-------------------------------------------\nERROR: " << strE
         << endl;
  }

  return 0;
}
