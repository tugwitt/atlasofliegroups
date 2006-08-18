/*!
\file
  \brief Definitions and declarations for the CartanClass and Fiber classes.
*/

/*
  This is cartanclass.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef CARTANCLASS_H  /* guard against multiple inclusions */
#define CARTANCLASS_H

#include "cartanclass_fwd.h"

#include "rootdata_fwd.h"
#include "tori_fwd.h"

#include "bitmap.h"
#include "gradings.h"
#include "latticetypes.h"
#include "partition.h"
#include "realform.h"
#include "setutils.h"
#include "size.h"
#include "subquotient.h"

namespace atlas {

/******** function declarations **********************************************/

namespace cartanclass {

  void compactTwoRho(latticetypes::Weight&, unsigned long, const Fiber&, 
		     const rootdata::RootDatum&);

  void restrictGrading(gradings::Grading&, const rootdata::RootSet&,
		       const rootdata::RootList&);

  void specialGrading(gradings::Grading&, const Fiber&, realform::RealForm, 
		      const rootdata::RootDatum&);

  void toMostSplit(rootdata::RootList&, const Fiber&, realform::RealForm,
		   const rootdata::RootDatum&);
}

/******** type definitions ***************************************************/

namespace cartanclass {

  /*!
  \brief Describes the fiber (over a fixed involution of a torus) in
  twisted Tits group.

  We fix an involutive automorphism tau of the torus H, and consider
  the collection of all strong real forms x of G that induce tau on H,
  modulo the conjugation action of H. Such a strong real form has
  square equal to some element z in Z(G)^delta.  Changing the element
  z by (1+delta)Z is innocuous, and there are finitely many
  possibilities for Z^delta/[(1+delta)Z].  For each fixed value of z,
  this collection of strong real forms (modulo H conjugation) has a
  simply transitive action (left multiplication) of the fiber group
  H^{-tau}/(1-tau)H.  The fiber group is a Z/2Z vector space.
  
  Equivalence classes of strong real forms (still with a fixed value
  of z) give a partition of the fiber group (depending on z).  These
  partitions are stored in d_strongReal, and accessed by the function
  strongReal.

  An equivalence class of strong real forms may be labelled by a pair
  of integers: the second labelling the element z (or rather its coset
  modulo [(1+delta)Z]), and the first the equivalence class in the
  fiber group.
  */
class Fiber {

 protected:

  /*!
  \brief Pointer to a torus defined over R.
  
  Represented as the lattice Z^n endowed with an involutive
  automorphism (represented by its n x n integer matrix).
  */
  tori::RealTorus* d_torus;

  /*!
  \brief RootSet flagging the complex roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the complex roots
  (those roots alpha with tau(alpha) equal neither to alpha nor to
  -alpha).
  */
  rootdata::RootSet d_complex;

  /*!
  \brief RootSet flagging the imaginary roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the imaginary
  roots (those roots alpha with tau(alpha)=alpha).
  */
  rootdata::RootSet d_imaginary;

  /*!
  \brief RootSet flagging the real roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the real roots
  (those roots alpha with tau(alpha)=-alpha).
  */
  rootdata::RootSet d_real;
  
  /*!
  \brief RootList holding the numbers of the simple imaginary
  roots.  

  These are simple for the positive imaginary roots given by the
  (based) RootDatum.  They need not be simple in the entire root
  system.
  */
  rootdata::RootList d_simpleImaginary;

  /*!
  \brief Records the action of the Cartan involution tau
  on the set of roots.

  root \#j is sent to root \#d_rootInvolution[j].
  */
  setutils::Permutation d_rootInvolution;
  
  /*!
  \brief Fiber group.  
  
  In terms of the complex torus H and the Cartan involution tau, it
  is equal to F=H^{-tau}/(1-tau)H: the group of fixed points of
  -tau, modulo its identity component.  

  Recall that the group of real points H(R) is a product of p circles,
  q copies of R^x, and r copies of C^x.  It is easy to see that
  F=(Z/2Z)^p, with one factor coming from each circle.   
  
  The fiber group is represented as a subquotient of the group H(2) of
  elements of order 2 in the torus.  Write Y for the lattice of one
  parameter subgroups of H; then H(2)=Y/2Y naturally.  The Cartan
  involution tau is always represented by the matrix q of its action on
  the character lattice X, which is dual to Y; so the action of tau on
  Y is given by the matrix q^t.  The action of -tau on Y is by -q^t.

  Because H is the product (not direct) of (H^tau)_0 and (H^{-tau})_0,
  every component of H^{-tau} meets H^tau.  The intersection of these
  two groups consists of H(2)^tau, the tau-fixed elements of order 2.  It
  follows that the component group F may be identified with
  H(2)^tau/{H(2) cap (1-tau)H}.  This is how F is stored as a
  subquotient of H(2)=Y/2Y=(Z/2Z)^n.  

  The larger space is computed as the kernel of the action of (1+tau)
  on H(2); that is, as the kernel of the reduction mod 2 of the matrix
  I+q^t.  The smaller space {H(2) cap (1-tau)H} consists of the
  elements of order 2 in the connected subtorus (H^{-tau})_0.  It is
  therefore computed as the image in H(2) of the -1 eigenlattice
  Y^{-tau} (the lattice of one parameter subgroups of (H^{-tau})_0).
  */
  latticetypes::ComponentSubquotient d_fiberGroup;

  /*!
  \brief Fiber group for the adjoint group of G.

  Writing H_ad for the complex torus in G_ad, and still writing tau
  for the Cartan involution, this is F_ad=H_ad^{-tau}/(1-tau)H_ad. The
  adjoint covering G-->G_ad defines a natural map F-->F_ad, but this
  map may be neither injective nor surjective.

  The adjoint fiber is computed from the action of tau on the adjoint
  lattice of one-parameter subgroups Y_ad just as the fiber group is
  computed. The lattice Y_ad has as basis the fundamental weights.  
  */
  latticetypes::ComponentSubquotient d_adjointFiberGroup;

  /*!
  \brief Matrix (over Z/2Z) of the map from the fiber
  group for G to the fiber group for G_ad.  

  Although the map is induced by Y/2Y-->Y_ad/2Y_ad, this matrix is
  _not_ just the reduction mod 2 of the matrix of Y-->Y_ad.  Each
  subquotient F and F_ad is equipped (by the row reduction algorithms
  in NormalSubquotient) with a distinguished basis, and what is
  recorded in d_toAdjoint is the matrix with respect to those bases.
  In particular, the size of the matrix is (dim F_ad) x (dim F). 
  */
  latticetypes::ComponentMap d_toAdjoint;
  
  /*!
  \brief Stabilizer of the grading in the adjoint fiber
   group.

  Each real form defines a grading of the imaginary root system,
  obtained by pairing with the simple imaginary roots. It should be
  the case that the grading entirely defines the real form, i.e., the
  corresponding W_i-orbit. However, this doesn't mean that the map
  from real form parameters to gradings has to be injective. The
  grading group contains the kernel of this map at the level of the
  adjoint fiber group.

  DV: I think this group must always be trivial.  The dual statement
  is that in the real points of a simply connected group, the
  component group of a torus is generated by the m_alpha for alpha
  real. Fokko includes in the code (cartanclass.cpp, after the
  calculation of d_gradingGroup) an "assert" that d_gradingGroup is
  trivial.
  */
  latticetypes::ComponentSubspace d_gradingGroup;
  
  /*!
  \brief Images in fiber group of the simple imaginary coroots.
  */
  latticetypes::ComponentList d_mAlpha;
  
  /*!
  \brief Images in adjoint fiber group of the simple imaginary coroots.
  */
  latticetypes::ComponentList d_adjointMAlpha;
  
  /*!
  \brief Grading with all simple imaginary roots noncompact.

  The grading is a RankFlags; that is, a BitSet<RANK_MAX>.  It flags
  the noncompact imaginary roots among the simple imaginary roots; so
  all the bits up to imaginaryRank are 1.
  */
  gradings::Grading d_baseGrading;
  
  /*!
  \brief Grading \#j flags the simple imaginary roots whose
  grading is changed by canonical basis vector \#j in the adjoint
  fiber group.
  */
  gradings::GradingList d_gradingShift;
  
  /*!
  \brief Flags the noncompact imaginary roots for the base grading
  among all the roots.
  */
  rootdata::RootSet d_baseNoncompact;
  
  /*!  
  \brief RootSet \#j flags the imaginary roots whose grading is
  changed by canonical basis vector \#j of the adjoint fiber group.
  */
  rootdata::RootSetList d_noncompactShift;
  
  /*!
  \brief Partition of the adjoint fiber group according to weak real
  forms.

  The imaginary Weyl group acts on the adjoint fiber group; the
  partition is by orbits of this action.  It is constructed by the
  function weakReal in the Helper class of the unnamed namespace of
  cartanclass.cpp.
  */
  partition::Partition d_weakReal;
  
  /*!
  \brief Partition of the weak real forms according to the
  corresponding classes in Z(G)^delta/[(1+delta)Z(G)].

  Constructed by the function realFormPartition in the Helper
  class of the unnamed namespace of cartanclass.cpp; details in that
  documentation.
  */
  partition::Partition d_realFormPartition;
  
  /*!  
  \brief Partitions of Fiber group cosets corresponding to the
  possible square classes in Z^delta/[(1+delta)Z].

  The Fiber group acts in a simply transitive way on strong real forms
  (inducing tau on H) with a fixed square in Z^delta.  The number of
  squares that occur (modulo (1+delta)Z) is equal to the number c of
  classes in the partition d_weakReal.  The collection of strong real
  forms is therefore a collection of c cosets of the fiber group F.
  Each of these c cosets is partitioned into W_i orbits; these orbits
  are described by the c partitions in d_strongReal.
  */
  std::vector<partition::Partition> d_strongReal;
  
  /*!
  \brief Representative strong real form for each weak real form.
  
  A StrongRealFormRep is a pair of unsigned long.  The second indexes
  the value of the square of the strong real form in
  Z^delta/[(1+delta)Z].  The first indexes a W_im orbit on the
  corresponding coset of the fiber group.
  */
  std::vector<StrongRealFormRep> d_strongRealFormReps;

 public:

// constructors and destructors

  Fiber():d_torus(0)
    {}

  Fiber(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  virtual ~Fiber();

// copy and assignment

  Fiber(const Fiber&);

  Fiber& operator= (const Fiber&);

  void swap(Fiber&);

// accessors

/*!
  \brief Fiber group for the adjoint group of G.

  Writing H_ad for the complex torus in G_ad, and still writing tau
  for the Cartan involution, this is F_ad=H_ad^{-tau}/(1-tau)H_ad. The
  adjoint covering G-->G_ad defines a natural map F-->F_ad, but this
  map may be neither injective nor surjective.

  The adjoint fiber is computed from the action of tau on the adjoint
  lattice of one-parameter subgroups Y_ad just as the fiber group is
  computed. The lattice Y_ad has as basis the fundamental weights.  
*/
  const latticetypes::ComponentSubquotient& adjointFiberGroup() const {
    return d_adjointFiberGroup;
  }

/*!
\brief Dimension of the adjoint fiber group as a Z/2Z vector space.
*/
  size_t adjointFiberRank() const {
    return d_adjointFiberGroup.dimension();
  }

/*!
\brief Cardinality of the adjoint fiber group.
*/
  unsigned long adjointFiberSize() const {
    return d_adjointFiberGroup.size();
  }

  void compactRootSet(rootdata::RootSet&, unsigned long) const;

/*!
  \brief RootSet flagging the complex roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the complex roots
  (those roots alpha with tau(alpha) equal neither to alpha nor to
  -alpha).
*/
  const rootdata::RootSet& complexRootSet() const {
    return d_complex;
  }

/*!
  \brief Fiber group.  
  
  In terms of the complex torus H and the Cartan involution tau, it
  is equal to F=H^{-tau}/(1-tau)H: the group of fixed points of
  -tau, modulo its identity component.  

  Recall that the group of real points H(R) is a product of p circles,
  q copies of R^x, and r copies of C^x.  It is easy to see that
  F=(Z/2Z)^p, with one factor coming from each circle.   
  
  The fiber group is represented as a subquotient of the group H(2) of
  elements of order 2 in the torus.  Write Y for the lattice of one
  parameter subgroups of H; then H(2)=Y/2Y naturally.  The Cartan
  involution tau is always represented by the matrix q of its action on
  the character lattice X, which is dual to Y; so the action of tau on
  Y is given by the matrix q^t.  The action of -tau on Y is by -q^t.

  Because H is the product (not direct) of (H^tau)_0 and (H^{-tau})_0,
  every component of H^{-tau} meets H^tau.  The intersection of these
  two groups consists of H(2)^tau, the tau-fixed elements of order 2.  It
  follows that the component group F may be identified with
  H(2)^tau/{H(2) cap (1-tau)H}.  This is how F is stored as a
  subquotient of H(2)=Y/2Y=(Z/2Z)^n.  

  The larger space is computed as the kernel of the action of (1+tau)
  on H(2); that is, as the kernel of the reduction mod 2 of the matrix
  I+q^t.  The smaller space {H(2) cap (1-tau)H} consists of the
  elements of order 2 in the connected subtorus (H^{-tau})_0.  It is
  therefore computed as the image in H(2) of the -1 eigenlattice
  Y^{-tau} (the lattice of one parameter subgroups of (H^{-tau})_0).
*/
  const latticetypes::ComponentSubquotient& fiberGroup() const {
    return d_fiberGroup;
  }

  /*!
  \brief Dimension of the fiber group as a Z/2Z vector space.
  */
  size_t fiberRank() const {
    return d_fiberGroup.dimension();
  }
  
  /*!
  \brief Cardinality of the fiber group: 2^dimension.
  */
  unsigned long fiberSize() const {
    return d_fiberGroup.size();
  }

  
void grading(gradings::Grading&, unsigned long) const;

  unsigned long gradingRep(const gradings::Grading&) const;

/*!
  \brief RootSet flagging the imaginary roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the imaginary
  roots (those roots alpha with tau(alpha)=alpha).
*/
  const rootdata::RootSet& imaginaryRootSet() const {
    return d_imaginary;
  }

  const latticetypes::LatticeMatrix& involution() const;


  void mAlpha(latticetypes::Component&, const rootdata::Root&) const;

  size_t minusRank() const;

  void noncompactRootSet(rootdata::RootSet&, unsigned long) const;

/*!
  \brief Number of weak real forms containing this Cartan.

  This is the number of orbits of the imaginary Weyl group on the
  adjoint fiber group.
*/
  size_t numRealForms() const {
    return d_weakReal.classCount();
  }

  size_t plusRank() const;

/*!
  \brief Partition of the weak real forms according to the
  corresponding classes in Z(G)/[(1+delta)Z(G)].

  A weak real form (always containing our fixed real torus) is an
  orbit of W_i on the adjoint fiber group.
*/
  const partition::Partition& realFormPartition() const {
    return d_realFormPartition;
  }

/*!
  \brief RootSet flagging the real roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the real roots
  (those roots alpha with tau(alpha)=-alpha).
*/
  const rootdata::RootSet& realRootSet() const {
    return d_real;
  }

/*!
\brief Action of the Cartan involution on root \#j.
*/
  rootdata::RootNbr rootInvolution(rootdata::RootNbr j) const {
    return d_rootInvolution[j];
  }

/*!
  \brief RootList holding the numbers of the simple imaginary
  roots.  

  These are simple for the positive imaginary roots given by the
  (based) RootDatum.  They need not be simple in the entire root
  system.
*/
  const rootdata::RootList& simpleImaginary() const {
    return d_simpleImaginary;
  }

/*!
 \brief Partitions of Fiber group cosets corresponding to the
  possible square classes in Z^delta/[(1+delta)Z].

  The Fiber group acts in a simply transitive way on strong real forms
  (inducing tau on H) with a fixed square in Z^delta.  The number of
  squares that occur (modulo (1+delta)Z) is equal to the number c of
  classes in the partition d_weakReal.  The collection of strong real
  forms is therefore a collection of c cosets of the fiber group F.
  Each of these c cosets is partitioned into W_i orbits; these orbits
  are described by the c partitions in d_strongReal.
*/
  const partition::Partition& strongReal(size_t j) const {
    return d_strongReal[j];
  }

/*!
  \brief Representative strong real form for real form \#rf.
  
  A StrongRealFormRep is a pair of unsigned long.  The second indexes
  the value of the square of the strong real form in
  Z^delta/[(1+delta)Z].  The first indexes a W_im orbit on the
  corresponding coset of the fiber group.
*/
  const StrongRealFormRep& strongRepresentative(size_t rf) const {
    return d_strongRealFormReps[rf];
  }

/*!

*/
  unsigned long toAdjoint(unsigned long) const;

/*!

*/
  unsigned long toWeakReal(unsigned long, size_t) const;

/*!
  \brief Real torus defined over R.
  
  Represented as the lattice Z^n endowed with an involutive
  automorphism (represented by its n x n integer matrix).
*/
  const tori::RealTorus& torus() const {
    return *d_torus;
  }

/*!
  \brief Partition of the weak real forms according to the
  corresponding classes in Z(G)^delta/[(1+delta)Z(G)].
*/
  const partition::Partition& weakReal() const {
    return d_weakReal;
  }

 };

   /*!
   \brief Represents a single stable conjugacy class of Cartan subgroups.

   Mathematically this means the complex torus in the complex group,
   together with an involutive automorphism of this torus (the Cartan
   involution).  (More precisely, it is a W-conjugacy class of
   involutive automorphisms.)  As the class is now used in the
   software, the Cartan involution must be in the inner class
   specified by ComplexReductiveGroup.  Most of the interesting
   information is contained in the two underlying Fiber classes
   d_fiber and d_dualFiber.  First of all, that is where the Cartan
   involution lives (since the involution is needed to define the
   fibers).  But the fiber classes also record for which real forms
   this stable conjugacy class of Cartan subgroups is defined; and,
   given the real form, which imaginary roots are compact and
   noncompact.
   */
class CartanClass {

 private:

  /*!
  \brief Class of the fiber group H^{-tau}/[(1-tau)H] for this
  Cartan.  
  
  Elements (very roughly) correspond to possible extensions
  of the real form tau from H to G.
  */
  Fiber d_fiber;
  
  /*!
  \brief Class of the fiber group for the dual Cartan.

  Elements of the dual fiber group are characters of the group of
  connected components of the real points of H.
  */
  Fiber d_dualFiber;

  /*!
  \brief Roots simple for the "complex factor" of W^tau.

  The subgroup W^tau of Weyl group elements commuting with the Cartan
  involution tau has two obvious (commuting) normal subgroups: the
  Weyl group W^R of the real (that is, fixed by -tau) roots, and the
  Weyl group W^iR of the imaginary (that is, fixed by tau) roots.  But
  this is not all of W^tau: W^tau is a semidirect product of (W^R x
  W^iR) with a group W^C, the first factor being normal.  Here is how
  to describe W^C.

  Write RC (standing for "complex roots") for the roots orthogonal to
  (a) the sum of positive real roots, and also (b) the sum of positive
  imaginary roots.  It turns out that RC as a root system is the
  direct sum of two isomorphic root systems RC_0 and RC_1 interchanged
  by tau.  (There is no canonical choice of this decomposition.)  The
  group W^tau includes W^C, the diagonal subgroup of W(RC_0) x
  W(RC_1).  The list d_simpleComplex is the numbers of the simple
  roots for (a choice of) RC_0.  That is, the Weyl group W^C is
  isomorphic to the Weyl group generated by the reflections
  corresponding to the numbers in d_simpleComplex.
  */
  rootdata::RootList d_simpleComplex;
  
  /*!
  \brief Size of the W-conjugacy class of tau.
  
  The number of distinct involutions defining the same stable
  conjugacy class of Cartan subgroups.
  */
  size::Size d_orbitSize;

public:

// constructors and destructors

  CartanClass() 
    {}

  CartanClass(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  ~CartanClass() 
    {}

// copy and assignment: defaults are ok for copy and assignment

  void swap(CartanClass&);

// accessors
  
  /*!
  \brief Class of the fiber group for the dual Cartan.

  Elements of the dual fiber group are characters of the group of
  connected components of the real points of H.
  */
  const Fiber& dualFiber() const {
    return d_dualFiber;
  }
  
  /*!
  \brief Class of the fiber group H^{-tau}/[(1-tau)H] for this
  Cartan.  
  
  Elements (very roughly) correspond to possible extensions
  of the real form tau from H to G.
  */
  const Fiber& fiber() const {
    return d_fiber;
  }


/*!
  \brief RootSet flagging the imaginary roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the imaginary
  roots (those roots alpha with tau(alpha)=alpha).
*/
  const rootdata::RootSet& imaginaryRootSet() const {
    return d_fiber.imaginaryRootSet();
  }


  /*!
  \brief Returns the matrix of the involution on the weight lattice of
  the Cartan subgroup. 
  */
  const latticetypes::LatticeMatrix& involution() const {
    return d_fiber.involution();
  }

  bool isMostSplit(unsigned long) const;

  /*!
  \brief Number of weak real forms for the dual group containing the
  dual Cartan.
  */
  unsigned long numDualRealForms() const {
    return d_dualFiber.numRealForms();
  }
  /*!
  \brief Number of weak real forms containing this Cartan.

  This is the number of orbits of the imaginary Weyl group on the
  adjoint fiber group.
  */
  unsigned long numRealForms() const {
    return d_fiber.numRealForms();
  }

  /*!
  \brief Number of possible squares of strong real forms mod
  (1+delta)Z.

  This is the number of classes in the partition of weak real forms
  according to Z^delta/[(1+delta)Z]. 
  */
  unsigned long numRealFormClasses() const {
    return d_fiber.realFormPartition().classCount();
  }

  /*! 
  \brief Size of the W-conjugacy class of tau.

  The number of distinct involutions defining the same stable conjugacy
  class of Cartan subgroups. 
   */
   unsigned long orbitSize() const {
    return d_orbitSize.toUlong();
  } 

   /*!
   \brief RootSet flagging the real roots.

   That is, a bitmap whose set bits are those corresponding to the
   numbers (within the list of roots in RootDatum) of the real roots
   (those roots alpha with tau(alpha)=-alpha).
   */
  const rootdata::RootSet& realRootSet() const {
    return d_fiber.realRootSet();
  }

  /*!
  \brief Action of the Cartan involution on root \#j.
  */
  rootdata::RootNbr rootInvolution(rootdata::RootNbr j) const {
    return d_fiber.rootInvolution(j);
  }

  /*!
  \brief RootList holding the numbers of the simple imaginary roots.

  These are simple for the positive imaginary roots given by the
  (based) RootDatum. They need not be simple in the entire root
  system.
  */
  const rootdata::RootList& simpleImaginary() const {
    return d_fiber.simpleImaginary();
  }

  /*!
  \brief Roots simple for the "complex factor" of W^tau.

  The subgroup W^tau of Weyl group elements commuting with the Cartan
  involution tau has two obvious (commuting) normal subgroups: the
  Weyl group W^R of the real (that is, fixed by -tau) roots, and the
  Weyl group W^iR of the imaginary (that is, fixed by tau) roots. But
  this is not all of W^: W^tau is a semidirect product of (W^R x W^iR)
  with a group W^C, the first factor being normal. Here is how to
  describe W^C.

  Write RC (standing for "complex roots") for the roots orthogonal to
  (a) the sum of positive real roots, and also (b) the sum of positive
  imaginary roots. It turns out that RC as a root system is the direct
  sum of two isomorphic root systems RC_0 and RC_1 interchanged by
  tau. (There is no canonical choice of this decomposition.) The group
  W^tau includes W^C, the diagonal subgroup of W(RC_0) x W(RC_1). The
  list d_simpleComplex is the numbers of the simple roots for (a
  choice of) RC_0. That is, the Weyl group W^C is isomorphic to the
  Weyl group generated by the reflections corresponding to the numbers
  in d_simpleComplex.
   */
  const rootdata::RootList& simpleComplex() const {
    return d_simpleComplex;
  }
  
  /*!
  \brief RootList holding the numbers of the simple real roots.

  These are simple for the positive real roots given by the (based)
  RootDatum. They need not be simple in the entire root system. 
  */
  const rootdata::RootList& simpleReal() const {
    return d_dualFiber.simpleImaginary();
  }

  /*!
  \brief Partitions of Fiber group cosets corresponding to the
  possible square classes in Z^delta/[(1+delta)Z].

  The Fiber group acts in a simply transitive way on strong real forms
  (inducing tau on H) with a fixed square in Z^delta. The number of
  squares that occur (modulo (1+delta)Z) is equal to the number c of
  classes in the partition d_weakReal. The collection of strong real
  forms is therefore a collection of c cosets of the fiber group
  F. Each of these c cosets is partitioned into W_i orbits; these
  orbits are described by the c partitions in d_strongReal.
  */
  const partition::Partition& strongReal(size_t j) const {
    return d_fiber.strongReal(j);
  }

  /*!
  \brief Returns the image of x in the adjoint fiber group.

  Precondition: x is a valid element in the fiber group. 
  */
  unsigned long toAdjoint(unsigned long x) const {
    return d_fiber.toAdjoint(x);
  }

  /*!
  \brief Returns the class number in the weak real form partition of
  the strong real form \#c in real form class rfc.

  The pair (c,rfc) is the software representation of an equivalence
  class of strong real forms (always assumed to induce tau on H). The
  integer rfc labels an element of Z^delta/[(1+delta)Z], thought of as
  a possible square value for strong real forms. The fiber group acts
  simply transitively on strong real forms with square equal to
  rfc. The integer c labels an orbit of W_i on this fiber group coset;
  this orbit is the equivalence class of strong real forms.
  */
  unsigned long toWeakReal(unsigned long c, size_t j) const {
    return d_fiber.toWeakReal(c,j);
  }

  /*!  
  \brief Partition of the weak real forms according to the
  corresponding classes in Z(G)^delta/[(1+delta)Z(G)].
  */
  const partition::Partition& weakReal() const {
    return d_fiber.weakReal();
  }

};

}

}

#endif
