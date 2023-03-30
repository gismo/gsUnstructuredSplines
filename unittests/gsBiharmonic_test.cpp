/** @file gsBiharmonic_test.cpp

    @brief Tests for Approx. C1 basis with solving biharmonic problem.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#include "gismo_unittest.h"

using namespace gismo;

void runBiharmonicTest ()
{
  gsInfo << "Test loaded successful\n";
}

SUITE(gsBiharmonic_test)
{
  TEST(Approx_C1_test)
  {
    runBiharmonicTest();
  }
}
