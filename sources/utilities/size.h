/*!
\file
  This is size.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef SIZE_H  /* guard against multiple inclusions */
#define SIZE_H

#include <cstring>
#include "constants.h"

/******** type declarations **************************************************/

namespace atlas {

namespace size {

  template<typename C> class SizeType;

  // the following is safe for RankMax<=128 (for C128 there are 255 factors 2)
  typedef signed char BaseType;
  typedef unsigned char UnsignedBaseType;
  typedef SizeType<BaseType> Size;

}


/******** constant declarations **********************************************/

namespace size {

  /*! \brief
  A template to indicate the (manually computed) ordinal (position on the list
  of primes) of the largest prime factor of a Weyl group of rank at most n.

  The ordinal is recorded in PrimesMax<n>::value, and used to make the
  SizeType class for storing Weyl group orders. This class needs instantiation
  only for |n| a possible value of |constants::RANK_MAX|; there powers of 2 in
  the order of the size of a word in terms of bits are most useful.
  */
  template<unsigned long n> class PrimesMax
  {
  public:
    static const unsigned long value; // no value supplied in generic case
  };


  // predefined values for likely instances of RANK_MAX
  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most 8.

  This is the fourth prime 7.
  */
  template<> class PrimesMax<8> {
  public:
    static const unsigned long value = 4ul;
  };

  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most 16.

  This is the seventh prime 17 (appearing only in the Weyl group S_17
  of type A16).
  */
  template<> class PrimesMax<16> {
  public:
    static const unsigned long value = 7ul;
  };

  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most 32.

  This is the eleventh prime 31.
  */

  template<> class PrimesMax<32> {
  public:
    static const unsigned long value = 11ul;
  };

  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most 64.

  This is the eighteenth prime 61.
  */

  template<> class PrimesMax<64> {
  public:
    static const unsigned long value = 18ul;
  };

  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most RANK_MAX.

  With RANK_MAX=16, this is 7 (for the seventh prime 17).
  */

  const size_t PRIMES_MAX = PrimesMax<constants::RANK_MAX>::value;

}  // |namespace size|

/******** function declarations **********************************************/

namespace size {

  unsigned long prime(size_t);

  Size factorial (unsigned long n);

}

/******** type definitions ***************************************************/

namespace size {

  /*!
  \brief Stores a positive integer as product of prime powers, using
  the first PRIMES_MAX primes.

  The exponent of the jth prime is |d_exp[j]|.  The reason for using
  this is that the software must occasionally deal with integers too
  large to be unsigned long; mostly these appear as cardinalities of
  Weyl groups.  The constant |PRIMES_MAX| is chosen so that the
  cardinality of any Weyl group of rank at most RANK_MAX can be
  represented as a |SizeType|.
  */
template<typename C> class SizeType
{
  C d_exp[PRIMES_MAX];

 public:
// constructors and destructors
  SizeType() {
    std::memset(d_exp,0,PRIMES_MAX*sizeof(C));
  }

  explicit SizeType(unsigned long);

  ~SizeType()
    {}

// copy and assignment

// (defaults suffice)

// accessors
  C operator[] (size_t j) const { return d_exp[j]; }

  bool operator== (const SizeType& c) const {
    return std::memcmp(d_exp,c.d_exp,PRIMES_MAX*sizeof(C))==0;
  }

  bool operator!= (const SizeType& c) const {
    return std::memcmp(d_exp,c.d_exp,PRIMES_MAX*sizeof(C))!=0;
  }

  bool hasOverflow() const;

  bool hasOverflow(size_t) const;

  unsigned long piece(size_t) const;

  unsigned long toUlong() const;

// manipulators
  C& operator[] (size_t j) {
    return d_exp[j];
  }

  SizeType& operator*= (const SizeType&);

  SizeType& operator*= (unsigned long n) { return operator*=(SizeType(n)); }

  SizeType& operator/= (const SizeType&);

  void reset() {
    std::memset(d_exp,0,PRIMES_MAX*sizeof(C));
  }

  void twoShift(C n) { // multiplication by 2^n; prime #0 is 2
    d_exp[0] += n;
  }
};  // |template<typename C> class SizeType|

}  // |namespace size|

}  // |namespace atlas|

#include "size_def.h"

#endif
