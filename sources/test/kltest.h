/*!
\file
\brief Declaration of the functions in namespace kltest.

These functions are designed to test a mathematical assertion used in
the implementation of the KL algorithm.
*/
/*
  This is kltest.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef KLTEST_H  /* guard against multiple inclusions */
#define KLTEST_H

#include "kgb_fwd.h"
#include "kl_fwd.h"

#include "permutations.h"

namespace atlas {

/******** type declarations *************************************************/

namespace kltest {

}

/******** function declarations **********************************************/

namespace kltest {

  bool checkBasePoint(const kgb::KGB&);

  void dualityPermutation(permutations::Permutation&, const kl::KLContext&);

  bool dualityVerify(const kl::KLContext& klc, const kl::KLContext& dual_klc);
}

/******** type definitions **************************************************/

namespace kltest {

}

}

#endif
