/*!
\file
\brief Class definition for the class DescentStatus.
*/

/*
  This is descents.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef DESCENTS_H  /* guard against multiple inclusions */
#define DESCENTS_H

#include <cstring>

namespace atlas {

/******** constant declarations *********************************************/

namespace descents {}

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace descents {


  /*!
\brief Describes the descent status of each simple root for a single
representation.

    There are eight possibilities for the descent status of a representation
    parameter w.r.t. a simple reflection, so this information could be packed
    in three bits. However, this would lead to having some packets lie across
    word boundaries, with the ensuing complications. Four bits is a good
    choice; here we take the lazy way of using even eight bits, as this makes
    the reading even a bit easier. We might come back to four at some later
    change---this should not require a change in user interface.
  */
class DescentStatus {

  /*!
\brief Value of byte \#j specifies the descent status of simple root
\#j.

Value should be 0 through 7; values 0 through 3 are ascents, and 4
through 7 are descents.
  */
  unsigned char d_data[constants::RANK_MAX];
  static const unsigned ValMask = constants::charBits - 1;

  /*!
\brief Bitwise "and" of Value with this is non-zero if Value is one of
ImaginaryCompact, ComplexDescent, RealTypeII, or RealTypeI (numbers 4--7)
  */
  static const unsigned DescentMask = 0x4ul;

  /*!
\brief Bitwise "and" of Value with this is equal to this if Value is either
ComplexDescent (5) or RealTypeI (7)
  */
  static const unsigned DirectRecursionMask = 0x5ul;

 public:

  enum Value { ComplexAscent, RealNonparity, ImaginaryTypeI, ImaginaryTypeII,
	       ImaginaryCompact, ComplexDescent, RealTypeII, RealTypeI };

  /*!
\brief Tests whether Value is 4 through 7.  These are the descents.

The simple roots passing this test comprise the tau invariant for the
representation.
  */
  static bool isDescent(Value v) {
    return (v & DescentMask)!=0;
  }

  /*!
\brief Tests whether both bits of DirectRecursionMask are set

In the case of a complex descent or a real type I descent there is a simple
recursion formula for the KL element.
  */
  static bool isDirectRecursion(Value v) {
    return (v & DirectRecursionMask) == DirectRecursionMask;
  }

// constructors and destructors
  DescentStatus() { // sets statuses of all simple roots to 0 (ComplexAscent)
    std::memset(d_data,0,constants::RANK_MAX);
  }

  ~DescentStatus() {}

// copy and assignment (these copy statuses of all simple roots)
  DescentStatus(const DescentStatus& ds) {
    std::memcpy(d_data,ds.d_data,constants::RANK_MAX);
  }

  DescentStatus& operator=(const DescentStatus& ds) {
    std::memcpy(d_data,ds.d_data,constants::RANK_MAX);
    return *this;
  }

// accessors

/*!
\brief Returns descent status of simple root \#s.
*/
  Value operator[] (size_t s) const {
    return static_cast<Value> (d_data[s]); // cast converts integer to enum
  }

// manipulators

/*!
\brief Sets the descent status of simple root \#s to v.
*/
  void set(size_t s, Value v) {
    d_data[s] = v; // no cast needed here; enum value converts to integral type
  }
};//class DescentStatus
  /*!
\brief Describes the hdescent status of each simple root (for the
small delta-fixed Hecke algebra) for a single h-fixed representation.
The simple root is a delta-orbit of simple roots for the big Hecke
algebra. If that simple root is delta-fixed, everything looks the same
as classically; I call those eight cases one-complex ascent, one-real nonparity,
one-noncompact imaginary type I, and so on. DON'T USE 8-15?

If alpha and delta(alpha) are orthogonal, a delta-fixed block element
gives theta commuting with delta, there are (perhaps?) 12
possibilities possibilities. Could be 1 or 2 or 4 discrete series?
Weyl group could be either -trivial (four discrete series, delta flips
the two middle ones?)  -s_alpha s_delta(alpha) (two discrete series,
both h-fixed) -four (one discrete series, h-fixed)

ADD SIXTEEN TO ALL THESE NUMBERS?

0) two-complex ascent (theta(alpha) positive, not in {\pm alpha, \pm
delta alpha}

1) two-semiimaginary ascent (theta(alpha)=delta(alpha); case of complex G,
involution commuting with s_alpha)

2) two-imaginary noncompact type I-I ascent (four discrete series, two
delta-fixed)

3) two-imaginary noncompact type II-II ascent (one discrete series)

4) two-imaginary noncompact type I-II ascent (two discrete series,
both delta-fixed)

5) two-real nonparity "ascent" (irreducible principal series on
product of two split A1s; not strict ascent)

8) two-complex descent

9) two-semireal descent (complex, length down one)

10) two-real type II-II descent (four fin diml, two delta-fixed)

11) two-real type I-I descent (one fin diml)

12) two-real type II-I descent (two delta-fixed fin diml)

13) two-imaginary compact "descent") (product of two compact A1s; not
strict descent)

If alpha and delta(alpha) are adjacent, then Lusztig's generator is
length three in the original Hecke algebra; so call these cases
"three-*". NUMBER THEM 32+..?

0) three-complex ascent (example is #3 in SL(5) big block; length
three higher)

1) three-semiimaginary ascent (first complex, then double-valued
Cayley; including other order gives three terms length two higher, one
delta-fixed. Example is #0 or #24 in SL(5) big block.)

2) three-real nonparity (both rn; example #0 in dual to SU(4,1))

8) three-complex descent (example is #31 in SL(5))

9) three semireal descent (first single-valued inverse Cayley, then C-)

10) three-compact imaginary (not possible when delta corrs to theta?)

For now I am treating only complex groups, so there are four
possibilities: 

0) two-complex ascent (theta(alpha) positive, not in {\pm alpha, \pm
delta alpha}

1) two-semiimaginary ascent (theta(alpha)=delta(alpha); case of complex G,
involution commuting with s_alpha)

2) two-complex descent

3) two-semireal descent (complex, length down one)

  */
class hDescentStatus {

  /*!
\brief Value of byte \#j specifies the h-descent status of h-simple root
\#j.

Value should be 0 through 3; values 0 through 1 are ascents, and 2
through 3 are descents.  (Eventually at least 26 values as described
above. I don't understand this constant, but it isn't used.)
  */
  unsigned char d_data[constants::RANK_MAX];
  static const unsigned hValMask = constants::charBits - 1;

  /*!
\brief Bitwise "and" of hValue with this is non-zero if hValue is one of
TwoComplexDescent or TwoSemiRealDescent (numbers 2-3)
  */
  static const unsigned hDescentMask = 0x2ul;

  /*!
\brief Bitwise "and" of hValue with this is equal to this if hValue is either
TwoComplexDescent (2) or TwoSemiRealDescent (3)
  */
  static const unsigned DirecthRecursionMask = 0x2ul;

 public:

  enum hValue { TwoComplexAscent, TwoSemiImaginary, TwoComplexDescent, 
	       TwoSemiReal };

  /*!
\brief Tests whether Value is 2 through 3.  These are the descents.

The simple roots passing this test comprise the tau invariant for the
representation.
  */
  static bool ishDescent(hValue v) {
    return (v & hDescentMask)!=0;
  }

  /*!
\brief Tests whether bit of DirecthRecursionMask is set

Maybe should separate the cases when division by q+1 is needed?
  */
  static bool isDirecthRecursion(hValue v) {
    return (v & DirecthRecursionMask) == DirecthRecursionMask;
  }

// constructors and destructors
  hDescentStatus() { // sets statuses of all simple roots to 0
		     // (TwoComplexAscent) 
    std::memset(d_data,0,constants::RANK_MAX);
  }

  ~hDescentStatus() {}

// copy and assignment (these copy statuses of all simple roots)
  hDescentStatus(const hDescentStatus& ds) {
    std::memcpy(d_data,ds.d_data,constants::RANK_MAX);
  }

  hDescentStatus& operator=(const hDescentStatus& ds) {
    std::memcpy(d_data,ds.d_data,constants::RANK_MAX);
    return *this;
  }

// accessors

/*!
\brief Returns descent status of simple root \#s.
*/
  hValue operator[] (size_t s) const {
    return static_cast<hValue> (d_data[s]); // cast converts integer to enum
  }

// manipulators

/*!
\brief Sets the descent status of simple root \#s to v.
*/
  void set(size_t s, hValue v) {
    d_data[s] = v; // no cast needed here; enum value converts to integral type
  }
};

}//namespace descents

}

#endif
