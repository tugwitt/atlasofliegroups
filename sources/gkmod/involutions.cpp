/*!
\file
\brief Implementation of the class InvolutionSet.
*/
/*
  This is involutions.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "involutions.h"

#include <map>
#include <set>
#include <stack>

#include "complexredgp.h"
#include "latticetypes.h"
#include "rootdata.h"

namespace atlas {

  // namespace {
  namespace involutions {
  namespace helper {

  using namespace involutions;

      /*!
\brief Derived class of InvolutionSet, to carry out the construction
of InvolutionSet.
      */
  class Helper:public InvolutionSet {

  private:

    /*!
\brief Entry s is the inner representation in the dual
  Weyl group of the generator with inner representation s.

They both have the same outer representation.
    */
    weyl::WeylInterface d_toDualWeyl;

  public:

    // constructors and destructors
    explicit Helper(complexredgp::ComplexReductiveGroup&);

    virtual ~Helper() {}

    // manipulators
    void fill(const complexredgp::ComplexReductiveGroup&);

    void fillCartan(const complexredgp::ComplexReductiveGroup&);

    void fillDualInvolutions(const weyl::WeylGroup&);

    void weylCorrelation(const complexredgp::ComplexReductiveGroup&);
 };

  }

  }


/*****************************************************************************

        Chapter I -- The InvolutionSet class

  ... explain here when it is stable ...

******************************************************************************/

namespace involutions {

  using namespace atlas::involutions::helper;

InvolutionSet::InvolutionSet()

{}

InvolutionSet::InvolutionSet(complexredgp::ComplexReductiveGroup& G)

/*!
\brief Constructor; computes the set of twisted involutions for G.

  NOTE: the set of involutions is actually constructed twice: once by
  fillCartan(), and then here. This is unfortunate, but negligible in the
  grand scheme of things, and having G remember its involution set runs
  against the principle of keeping the ComplexReductiveGroup class fairly
  small. Also this should be easy enough to change in the future if desired.
*/

{
  Helper help(G);
  swap(help);
}

/******** copy, assignment and swap *****************************************/
void InvolutionSet::swap(InvolutionSet& other)

{
  std::swap(d_size,other.d_size);
  std::swap(d_rank,other.d_rank);
  d_action.swap(other.d_action);
  d_cartan.swap(other.d_cartan);
  d_involution.swap(other.d_involution);
  d_dualInvolution.swap(other.d_dualInvolution);

  return;
}

/******** assignment *******************************************************/
size_t InvolutionSet::involutionNbr(const weyl::WeylElt& w, 
				    const weyl::WeylGroup& W) const

/*!
\brief Index of w in d_involution.

  Precondition: w is a twisted involution for W;

  Algorithm: find a reduced expression of w as an involution, and follow the
  action pointers.
*/

{
  using namespace weyl;

  WeylWord ww;
  W.involutionOut(ww,w);

  size_t x = 0;

  for (size_t j = 0; j < ww.size(); ++j)
    x = action(ww[j],x);

  return x;
}

}

/*****************************************************************************

        Chapter II -- The Helper class

  ... explain here when it is stable ...

******************************************************************************/

// namespace {
  namespace involutions {
  namespace helper {

Helper::Helper(complexredgp::ComplexReductiveGroup& G)

{
  using namespace weyl;

  // fill in the full cartan structure for G
  G.fillCartan();

  d_size = G.numInvolutions();
  d_rank = G.semisimpleRank();

  weylCorrelation(G);
  fill(G);
}

/******** manipulators *******************************************************/

void Helper::fill(const complexredgp::ComplexReductiveGroup& G)

/*!
\brief Fills the tables.

  Precondition: size and rank have been set; W is the relevant Weyl group.

  Explanation: action(s,w) is s.w if s and w twisted-commute, s.w.twist(s) 
  otherwise.

  Algorithm: we fill the table in order of increasing involution length, by the
  naive algorithm of looking at all the slots which have not yet been filled,
  computing the result and looking it up in a set of elements for the next
  length.
*/

{  
  using namespace weyl;

  const WeylGroup& W = G.weylGroup();

  // fill the action and involution tables
  d_action.resize(d_rank);
  for (size_t s = 0; s < d_action.size(); ++s)
    d_action[s].resize(d_size,UndefInvolution);

  d_involution.resize(d_size);

  std::map<WeylElt,size_t> found;
  size_t nextLength = 1;
  size_t firstNew = 1;

  for (size_t x = 0; x < d_size; ++x) {

    if (x == nextLength) { // update
      nextLength += found.size();
      found.clear();
    }

    for (size_t s = 0; s < d_rank; ++s) {

      if (action(s,x) != UndefInvolution)
	continue;

      const WeylElt& w = involution(x);
      WeylElt sw = w;

      if (W.hasTwistedCommutation(s,w)) { // action is product
	W.leftProd(sw,s);
      } else { // action is twisted commutation
	W.twistedConjugate(sw,s);
      }

      if (found.insert(std::make_pair(sw,firstNew)).second) { 
	// found a new element
	d_involution[firstNew] = sw;
	d_action[s][x] = firstNew;
	d_action[s][firstNew] = x;
	++firstNew;
      } else { // sw is already known
	size_t sx = found.find(sw)->second;
	d_action[s][x] = sx;
	d_action[s][sx] = x;
      }
    }
  }

  // fill in the dual involution table
  fillDualInvolutions(W);

  // fill in the cartan table
  fillCartan(G);

  return;
}

void Helper::fillCartan(const complexredgp::ComplexReductiveGroup& G)

/*!
\brief Fills the Cartan table.

  Precondition: the action and involution tables have been filled;

  Explanation: it is important that we use the same numbering of Cartan
  subgroups as in G. The algorithm is to use the representative of Cartan
  #j returned by G.twistedInvolution(j), locate that in the involution set,
  and then number its cross-orbit with j's.
*/

{
  using namespace weyl;

  const WeylGroup& W = G.weylGroup();

  d_cartan.resize(d_size);

  std::set<size_t> found;
  std::stack<size_t> toDo;

  for (size_t cn = 0; cn < G.numCartanClasses(); ++cn) {

    const WeylElt& w = G.twistedInvolution(cn);
    size_t x0 = involutionNbr(w,W);

    // find cross-orbit of x0
    found.clear();

    found.insert(x0);
    toDo.push(x0);

    while (!toDo.empty()) {
      
      size_t x = toDo.top();
      toDo.pop();
      
      for (Generator s = 0; s < d_rank; ++s) {
	const WeylElt& w = involution(x);
	if (W.hasTwistedCommutation(s,w))
	  continue;
	size_t sx = action(s,x);
	if (found.insert(sx).second)
	  toDo.push(sx);
      }
    }

    // write result
    std::set<size_t>::iterator found_end = found.end();
    for (std::set<size_t>::iterator i = found.begin(); i != found_end; ++i)
      d_cartan[*i] = cn;

  }

  return;
}

void Helper::fillDualInvolutions(const weyl::WeylGroup& W)

/*!
\brief Fills the dual involution table.

  Precondition: the action and involution tables have been filled; d_toDualWeyl
  is set;
*/

{
  using namespace weyl;

  d_dualInvolution.resize(d_size);

  for (size_t j = 0; j < d_size; ++j) {
    WeylElt v = W.longest();
    WeylElt w = involution(j);
    W.translate(w,d_toDualWeyl);
    W.twist(w);
    W.prod(v,w);
    W.invert(v);
    d_dualInvolution[j] = v;
  }
  
  return;
}

void Helper::weylCorrelation(const complexredgp::ComplexReductiveGroup& G)

/*!
\brief Fills in d_toDualWeyl.

  Explanation: d_toDualWeyl[s] is the inner representation in the dual
  Weyl group of the generator with inner representation s (they both
  have the same outer representation.)
*/

{  
  using namespace latticetypes;
  using namespace rootdata;
  using namespace weyl;

  const WeylGroup& W = G.weylGroup();

  // make dual Weyl group
  LatticeMatrix c;
  cartanMatrix(c,G.rootDatum());
  c.transpose();
  Twist t;
  W.outerTwist(t);
  WeylGroup dW(c,&t);

  // fill in d_toDualWeyl
  for (size_t s = 0; s < W.rank(); ++s) {
    WeylElt w;
    W.prod(w,s);
    WeylWord ww;
    dW.out(ww,w);
    d_toDualWeyl[s] = ww[0];
  }

  return;
}

}

}

}



