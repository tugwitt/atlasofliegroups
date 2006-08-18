/*!
\file 
\brief Class definitions and function declarations for the classes
TitsGroup and TitsElt.
*/
/*
  This is tits.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef TITS_H  /* guard against multiple inclusions */
#define TITS_H

#include "tits_fwd.h"

#include "rootdata_fwd.h"

#include "bitvector.h"
#include "constants.h"
#include "latticetypes.h"
#include "weyl.h"

namespace atlas {

/******** type declarations *************************************************/

namespace tits {

  /*!
\brief Element of (Z/2Z)^n, representing an element of T(2).

  The cocharacter lattice of T is always identified with Z^n; applying a
  cocharacter to -1 identifies the group T(2) of elements of order 2
  with (Z/2Z)^n.
  */ 
  typedef bitvector::BitVector<constants::RANK_MAX> TorusPart;

  /*!
\brief Square matrix of elements of Z/2Z, representing an endomorphism of T(2).
  */ 
  typedef bitvector::BitMatrix<constants::RANK_MAX> TorusMatrix;

}

/******** function declarations *********************************************/

namespace tits {

  void makeTwist(weyl::Twist&, const latticetypes::LatticeMatrix&,
		 const rootdata::RootDatum&);

}

/******** type definitions **************************************************/

namespace tits {

  /*!
\brief Represents an element of a Tits group.

An element is always written w.t, with w the canonical representative
in the Tits group of a Weyl group element w and t in T(2).
  */
class TitsElt {

 private:
  /*!
\brief Canonical Weyl part of the Tits group element.

If the Weyl group element w has a reduced decomposition s_1...s_r,
then the canonical representative is the product of the corresponding
Tits group elements sigma_1...sigma_r.  Because the sigma_i satisfy
the braid relations, this canonical representative is independent of
the choice of reduced decomposition. 
  */
  weyl::WeylElt d_w;

  /*!
\brief Factor in T(2) of the Tits group element.
  */
  TorusPart d_t;
  
 public:

// constructors and destructors

/*! 
\brief Constructs the identity element in the group
*/
  explicit TitsElt(size_t n):d_t(n) {} 

  /*!
\brief Constructs the Tits element with Weyl part w and torus factor
n.

The unsigned long n is regarded (in binary representation) as an
element of (Z/2Z)^RANK_MAX, which contains T(2).  Making RANK_MAX
exceed the machine word size would require a different interface for
this constructor. [DV 7/23/06.]
  */
  TitsElt(const weyl::WeylElt& w, size_t n):d_w(w),d_t(n) {}

// copy and assignment can be left to their defaults

// accessors
  bool operator< (const TitsElt& a) const {
    if (d_w == a.d_w)
      return d_t < a.d_t;
    else
      return d_w < a.d_w;
  }

  /*!
\brief Factor in T(2) of the Tits group element.
  */
  const TorusPart& t() const {
    return d_t;
  }

  /*!
\brief Canonical Weyl part of the Tits group element.
  */
  const weyl::WeylElt& w() const {
    return d_w;
  }

// manipulators

  /*!
\brief Multiplies the Tits group element on the right by x in T(2).
  */
  TitsElt& operator+= (const TorusPart& x) {
    d_t += x;
    return *this;
  }

  /*!
\brief Factor in T(2) of the Tits group element.
  */
  TorusPart& t() {
    return d_t;
  }

  /*!
\brief Canonical Weyl part of the Tits group element.
  */
  weyl::WeylElt& w() {
    return d_w;
  }

};

/*!
\brief Represents a finite subgroup of the normalizer in G of the Cartan T.
  
  We use a slight variant of the Tits group (also called extended Weyl
  group) as defined in J. Tits, J. of Algebra 4 (1966), pp. 96-116.

  The slight variant is that we include all elements of order two in the torus,
  instead of just the subgroup generated by the m_alpha (denoted h_alpha in
  Tits's paper.) Tits' original group may be defined by generators 
  sigma_alpha, alpha simple, subject to the braid relations and sigma_alpha^2=
  m_alpha; to get our group we just add a basis of elements of T(2) as 
  additional generators, and express the m_alpha in this basis.

  On a practical level, because the sigma_alpha satisfy the braid relations,
  any element of the Weyl group has a canonical lifting in the Tits group;
  so we may uniquely represent any element of the Tits group as a pair (w,t),
  with t in T(2) and w in W. The multiplication rules have to be thoroughly
  changed, though, to take into account the new relations.

  We have not tried to be optimally efficient here, as it is not expected that
  Tits computations will be significant computationally.

  Note on independence of choices: given a root alpha of T in G, the
  corresponding homomorphism phi_alpha from SL(2) to G is defined only
  up to conjugation by T.  This means that the generator sigma_alpha of
  the Tits group appears to be defined only up to multiplication by
  the image of the coroot alpha^vee.  But we are fixing a pinning,
  which means in particular that the maps phi_alpha for alpha simple
  are canonically defined (by the requirement that the standard
  pinning of SL(2) be carried into the pinning for G).  This means
  that the generator sigma_alpha (still for alpha simple) is canonically
  defined.
*/
class TitsGroup {

 private:

  /*!
\brief Pointer to the underlying Weyl group.
  */
  weyl::WeylGroup* d_weyl; // the underlying Weyl group
  
  /*!
\brief Dimension of the Cartan T.
  */
  size_t d_rank;
  
  /*!
\brief List of the images in character lattice mod 2 of the simple
roots.

Regarded as elements of order two in the dual torus ^vee T, these are
the elements m_alpha^vee.
  */
std::vector<TorusPart> d_simpleRoot;

  /*!
\brief List of the elements m_alpha (for alpha simple) in T(2).
  */
  std::vector<TorusPart> d_simpleCoroot;

  /*!
\brief Reduction mod 2 of the matrix of the defining involution of the
inner class.

Gives the action of this involution delta on T(2), for computing in
the delta coset of the Tits group.
  */
  TorusMatrix d_involution;

  /*!
\brief Permutation of the simple generators of the Tits group given by
the involution delta.

This is an array of unsigned char (labelling the generators) of size
RANK_MAX, giving the permutation of the generators induced by delta.
  */
  weyl::Twist d_twist;

// copy and assignment
// reserve and implement when necessary
  TitsGroup(const TitsGroup&);
  TitsGroup& operator= (const TitsGroup&);

 public:

// constructors and destructors
  TitsGroup() {}

  TitsGroup(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  ~TitsGroup();

// accessors
/*!
\brief Conjugate the TitsElt a by the generator for simple root \#s.

Note that the inverse of the generator sigma_alpha is sigma_alpha.m_alpha.
*/
  void conjugate(TitsElt& a, weyl::Generator s) const {

    leftProd(a,s);
    prod(a,s);
    a += d_simpleCoroot[s];
  }

  void invConjugate(TorusPart&, const weyl::WeylElt&) const;

  void leftProd(TitsElt&, weyl::Generator) const;

  /*!
\brief Length of the underlying Weyl group element.
  */
  unsigned long length(const TitsElt& a) const {
    return d_weyl->length(a.w());
  }

  void prod(TitsElt&, weyl::Generator) const;

  void prod(TitsElt&, const TitsElt&) const;

  /*!
\brief Rank of the underlying group.
  */
  const size_t rank() const {
    return d_rank;
  }

  /*!
\brief Applies to the element x of T(2) simple reflection s.

The reflection must be given its _outer_numbering (the Bourbaki one),
not the internal renumbering used by weylGroup.
  */
  void reflection(TorusPart& x, weyl::Generator s) const {
/*
  note: s must be an _outer_ generator
*/
    if (bitvector::scalarProduct(x,d_simpleRoot[s]))
      x += d_simpleCoroot[s];
  }

  /*!
\brief Element m_alpha of T(2) for simple coroot \#j.
  */
  const TorusPart& simpleCoroot(size_t j) const {
    return d_simpleCoroot[j];
  }

  /*!
\brief Image in the character lattice mod 2 of simple root \#j.
  */
  const TorusPart& simpleRoot(size_t j) const {
    return d_simpleRoot[j];
  }

  /*!
\brief Action of the defining involution of the inner class on simple
generator \#j.
  */
  size_t twist(size_t j) const {
    return d_twist[j];
  }

  /*!
\brief Twisted conjugates the TitsElt a by the generator for simple
root \#s.

This corresponds to conjugation of the coset a.delta, with delta the
defining involution of the inner class.  Note that the inverse of the
generator sigma_alpha is sigma_alpha.m_alpha.
  */
  void twistedConjugate(TitsElt& a, weyl::Generator s) const {
/*
  note: in the Tits group s^{-1} is s.m_s!
*/
    leftProd(a,s);
    prod(a,d_twist[s]);
    a += d_simpleCoroot[d_twist[s]];
  }

  const weyl::WeylGroup& weylGroup() const {
    return *d_weyl;
  }

// manipulators
  void swap(TitsGroup&);
};

}

}

#endif
