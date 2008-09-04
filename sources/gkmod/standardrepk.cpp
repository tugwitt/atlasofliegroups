/*!
\file
\brief Implementation for the classes
StandardRepK and KHatComputations.
*/
/*
  This is standardrepk.cpp

  Copyright (C) 2006, 2007 Alfred Noel
  With comments (C) 2008 by Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups version 0.2.4

  See file main.cpp for full copyright notice
*/

#include "standardrepk.h"

#include "cartanclass.h"
#include "complexredgp.h"
#include "realredgp.h"
#include "kgb.h"
#include "descents.h"
#include "lattice.h"
#include "rootdata.h"
#include "smithnormal.h"
#include "cartanset.h"
#include "tori.h"
#include "basic_io.h"
#include "intutils.h"
#include "tags.h"
#include "graph.h"
#include "prettyprint.h"
#include <iostream>
#include <stdexcept>

namespace atlas {



/*****************************************************************************

        Chapter I -- The StandardRepK class

******************************************************************************/

  namespace standardrepk {


// accessors

bool StandardRepK::operator< (const StandardRepK& rhs) const
{
  if (d_cartan != rhs.d_cartan) return d_cartan < rhs.d_cartan ;
  if (d_fiberElt != rhs.d_fiberElt) return d_fiberElt < rhs.d_fiberElt;
  return d_lambda < rhs.d_lambda;
}

bool StandardRepK::operator== (const StandardRepK& other) const
{
  return d_cartan == other.d_cartan
    and d_fiberElt == other.d_fiberElt
    and d_lambda == other.d_lambda;
  // when these components match, the restrictions to $K$ will be the same
}

} // namespace standardrepk


/*****************************************************************************

        Chapter II -- The KHatComputations class

******************************************************************************/

namespace standardrepk {

  // for now we only construct a KHatComputations from a KGB structure

KHatComputations::KHatComputations
  (const realredgp::RealReductiveGroup &GR, const kgb::KGB& kgb)
  : d_G(&GR.complexGroup())
  , d_KGB(kgb)
  , d_kHat()
  , d_standardHechtSchmid()
  , d_realForm(GR.realForm())
  , d_Tg(kgb::EnrichedTitsGroup::for_square_class(GR))
  , d_data(d_G->numCartanClasses())
{
  // the following loop should be restricted to Cartan classes for real form
  for (size_t r=0; r<d_data.size(); ++r)
  {
    // d_G->cartan[r] is rth CartanClass
    // which by now records the canonical involution for this class
    Cartan_info& ci=d_data[r];

    latticetypes::LatticeMatrix theta = d_G->cartan(r).involution();

    size_t n = rootDatum().rank();

    // put in $q$ the matrix of $\theta-1$
    latticetypes::LatticeMatrix q=theta;
    for (size_t j = 0; j < n; ++j)
      q(j,j) -= 1;

    // find Smith basis relative to $q$
    latticetypes::WeightList bs; matrix::initBasis(bs,n);
    latticetypes::CoeffList invf; smithnormal::smithNormal(invf,bs.begin(),q);

    latticetypes::LatticeMatrix basis(bs);
    latticetypes::LatticeMatrix inv_basis=basis.inverse();


    size_t l = invf.size();
    size_t f=0; while (f<l and invf[f]==1) ++f; // ignore invariant factors 1

    ci.torsionProjector.resize(l-f,n);

    for (size_t i=f; i<l; ++i)
    {
      assert(invf[i]==2); // other invariant factors must be 2

      using latticetypes::operator*=;
      ci.torsionLift.push_back(bs[i]);
      ci.torsionLift.back() *= 2; // double it, due to $(1/2)X^*$ basis

      for (size_t j=0; j<n; ++j)
	ci.torsionProjector.set_mod2(i,j,inv_basis(i,j));
    }
    if (l==n) // avoid construction from empty list
      ci.freeLift.resize(n,0); // so set matrix to empty rectangle
    else
      ci.freeLift= // or take the final |n-l| basis vectors as columns
	latticetypes::LatticeMatrix(&bs[l],&bs[n],tags::IteratorTag());

    ci.freeProjector.resize(n-l,n);

    for (size_t i =l; i<n; ++i) // copy final rows from inverse
      ci.freeProjector.copyRow(inv_basis,i-l,i); // row |i| to |i-l|


    latticetypes::LatticeElt e=inv_basis.apply(rootDatum().twoRho());
    for (size_t i=0; i<n; ++i)
      if (invf[i]!=2) e[i]=0; // clear components not on torsion part
    ci.p2rho=basis.apply(e);

    ci.fiber_modulus=latticetypes::SmallSubspace
      (latticetypes::SmallBitVectorList(tori::minusBasis(theta.transposed())),
       n);
  }
}

/******** accessors *******************************************************/

HCParam KHatComputations::project
  (size_t cn, const latticetypes::Weight& lambda) const
{
  using latticetypes::operator-=;
  using latticetypes::operator/=;
  const Cartan_info& ci=d_data[cn];

  latticetypes::Weight lambda_shifted=lambda;
  (lambda_shifted -= rootDatum().twoRho()) /= 2; // final division is exact
  return std::make_pair
    (ci.freeProjector.apply(lambda),
     ci.torsionProjector.apply(latticetypes::SmallBitVector(lambda_shifted))
     .data());
}

latticetypes::Weight KHatComputations::lift(size_t cn, HCParam p) const
{
  using latticetypes::operator+=;
  const Cartan_info& ci=d_data[cn];
  latticetypes::Weight result=ci.freeLift.apply(p.first); // lift free part

  result+=ci.p2rho; // offset to lift from torsion part, may have odd parts
  latticetypes::WeightList torsion_lift=ci.torsionLift;
  for (size_t i=0; i<torsion_lift.size(); ++i)
    if (p.second[i])
      result += torsion_lift[i]; // add even vectors representing torsion part

  return result;
}

StandardRepK KHatComputations::std_rep
  (const latticetypes::Weight& two_lambda, tits::TitsElt a) const
{
  weyl::TwistedInvolution sigma=a.tw();
  weyl::WeylElt w = d_G->cartanClasses().canonicalize(sigma);
  // now |sigma| is canonical and |w| conjugates |sigma| to |a.tw()|

  const weyl::WeylGroup& W=d_G->weylGroup();

  size_t cn = d_G->cartanClasses().classNumber(sigma);

  // conjugate towards canonical element
  {
    weyl::WeylWord ww=W.word(w);
    for (size_t i=0; i<ww.size(); ++i) // left-to-right for inverse conjugation
      d_Tg.basedTwistedConjugate(a,ww[i]);
  }
  assert(a.tw()==sigma); // we should now be at canonical twisted involution

  const rootdata::RootDatum& rd=rootDatum();
  latticetypes::Weight mu=W.imageByInverse(rd,w,two_lambda);
  StandardRepK result(cn,
		      d_data[cn].fiber_modulus.mod_image
		        (d_Tg.titsGroup().left_torus_part(a)),
		      project(cn,mu));

  const cartanclass::Fiber& cc=d_G->cartan(cn).fiber();
  for (size_t i=0; i<cc.imaginaryRank(); ++i)
    if (rd.scalarProduct(mu,cc.simpleImaginary(i))<0)
      return result; // without setting isStandard

  result.d_status.set(StandardRepK::IsStandard);
  return result;
}

StandardRepK KHatComputations::KGB_elt_rep(kgb::KGBElt z) const
{
  return std_rep(rootDatum().twoRho(),d_KGB.titsElt(z));
}

void KHatComputations::normalize(StandardRepK& stdrep) const
{
  if (not stdrep.isStandard())
  {
    std::cout << "cannot normalize properly continued standard representation"
	      << std::endl;
    return ;
  }

  rootdata::RootDatum rd = rootDatum();
  size_t cn=stdrep.d_cartan;
  const weyl::TwistedInvolution sigma_0 =
    d_G->cartanClasses().twistedInvolution(cn);

  const cartanclass::CartanClass& cc=d_G->cartan(cn);

  // lift lambda to X^*

  latticetypes::Weight lifted_lambda=lift(cn,stdrep.d_lambda);

  bitset::RankFlags bi_ortho_simples;
  { // find simple roots orthogonal to |real2rho| and |imaginary2rho|
    latticetypes::Weight real2rho=rd.twoRho(cc.realRootSet());
    latticetypes::Weight imaginary2rho=rd.twoRho(cc.imaginaryRootSet());
    for (size_t i=0; i<rd.semisimpleRank(); ++i)
      if (rd.isOrthogonal(real2rho,i) and rd.isOrthogonal(imaginary2rho,i))
	bi_ortho_simples.set(i);
  }

  const latticetypes::LatticeMatrix& theta = cc.involution();
  atlas::rootdata::RootSet cplx = cc.fiber().complexRootSet();

  //  go through the orthogonal list
  //  select the complex roots in the list

  bitset::RankFlags::iterator i;
  do
  {
    using latticetypes::operator+=;

    latticetypes::Weight tlifted_lambda=theta.apply(lifted_lambda);
    tlifted_lambda +=lifted_lambda;

    for (i= bi_ortho_simples.begin(); i(); ++i )
    {
      rootdata::RootNbr alpha = rd.simpleRootNbr(*i);
      rootdata::RootNbr talpha= cc.involution_image_of_root(alpha);

      if (cplx.isMember(alpha) and rd.scalarProduct(tlifted_lambda,alpha)<0)
      {
	assert (rd.isOrthogonal(alpha,talpha));
	rd.reflect(lifted_lambda,alpha);
	rd.reflect(lifted_lambda,talpha);
	break; // and continue do-while loop
      }
    } // for
  }
  while (i()); // while |for|-loop was exited through |break|


// put the new weakly dominant lambda back in stdrep

  stdrep.d_lambda=project(cn,lifted_lambda);
  stdrep.d_status.set(StandardRepK::IsNormal);

} // normalize


HechtSchmid
KHatComputations::HS_id(StandardRepK sr,rootdata::RootNbr alpha) const
{
  HechtSchmid id(sr);
  tits::TitsElt a=titsElt(sr);
  latticetypes::Weight lambda=lift(sr);
  latticetypes::Weight mu=rootDatum().reflection(lambda,alpha);

  if (d_Tg.grading(a,alpha))
  { // |alpha| is a non-compact root
    d_Tg.basedTwistedConjugate(a,alpha); // adds $m_\alpha$ to torus part
    id.add_lh(std_rep(mu, a));

    /* Now are equivalent:
       = |id.lh2->d_fiberElt==sr.d_fiberElt|
       = $m_s$ absorbed into quotient (i.e., $m_\alpha\in(X_*^- mod 2)$)
       = $\alpha^\vee( (X^*)^\theta ) = 2\Z$
       = $\alpha\notin(1-\theta_\alpha)X^*$  ($\theta_\alpha=\theta.s_\alpha$)
       = Cayley transform is type II
    */

    d_Tg.Cayley_transform(a,alpha);
    id.add_rh(std_rep(lambda, a));
    if (id.lh2->d_fiberElt==sr.d_fiberElt) // type II
    {
      using latticetypes::operator+=;
      latticetypes::Weight ra=rootDatum().simpleRoot(alpha);
      lambda += ra; // other possibility
      id.add_rh(std_rep(lambda, a));
      assert(id.rh1->d_lambda != id.rh2->d_lambda);
    }
  }
  else // |alpha| is a compact root; "easy" Hecht-Schmid identity
    // based twisted conjugation fixed the Tits element; just reflect weight
    id.add_lh(std_rep(mu, a)); // and no RHS

  return id;
}

atlas::latticetypes::LatticeCoeff
KHatComputations::product_simpleroot
  (const StandardRepK& s, rootdata::RootNbr k) const
{
  using latticetypes::operator+=;

  if (not s.isNormal())
    throw std::runtime_error("simpleroot: unnormalized standard rep");

  latticetypes::Weight lifted=theta_lift(s.d_cartan,s.d_lambda);

  return rootDatum().scalarProduct(lifted,k);
}

atlas::latticetypes::LatticeCoeff
KHatComputations::product_sumposroots(const StandardRepK& s) const
{
  using latticetypes::operator+=;

  if (not s.isNormal())
    throw std::runtime_error("sumposroots: unnormalized standard rep");

  latticetypes::Weight lifted=theta_lift(s.d_cartan,s.d_lambda);

  latticetypes::Weight dual2rho=rootDatum().dual_twoRho();

  return intutils::abs(latticetypes::scalarProduct(lifted,dual2rho));
}

PSalgebra
KHatComputations::theta_stable_parabolic
  (weyl::WeylWord& conjugator,
   const cartanset::CartanClassSet& cs,
   const size_t Cartan_nr) const
{
  const rootdata::RootDatum rd=cs.rootDatum();
  const weyl::WeylGroup W=cs.weylGroup();
  weyl::TwistedInvolution twi=cs.twistedInvolution(Cartan_nr);
  // latticetypes::LatticeMatrix theta=cs.involutionMatrix(twi);


/* Conjugate |twi| by simple complex roots to make the positive root system
   more theta stable. There is no hope for real roots, but for complex roots
   we try to achieve that the $\theta$-image of simple roots is positive. As
   our notion of positive roots is fixed, we conjugate $\theta$ itself (in
   fact twisted-conjugating |twi|) which changes the notions of
   imaginary/real/complex roots.
*/
  conjugator.resize(0);

  {
    size_t i;
    do // get list of simple complex roots
    {
      cartanclass::InvolutionData id(rd,cs.involutionMatrix(twi));
      for (i=0; i<rd.semisimpleRank(); ++i)
      {
	rootdata::RootNbr alpha=rd.simpleRootNbr(i);
	if (id.complex_roots().isMember(alpha)
	    and not rd.isPosRoot(id.root_involution()[alpha]))
	{
	  W.twistedConjugate(twi,i);
	  conjugator.push_back(i);
	  break; // and repeat do-while loop
	}
      } // for i
    }
    while (i!=rd.semisimpleRank());
  }
/*
   We have achieved that any real positive root is a sum of real simple roots.
   For consider a counterexample that is minimal: a real positive root
   $\alpha$ from which no real simple root can be subtracted while remaining
   positive. Then this must be a sum of simple roots that are either complex
   or imaginary; the $\theta$-images of such simple roots are all positive,
   but their sum is $-\alpha$ which is negative, a contradiction.

   Meanwhile |conjugator| is such that |W.twistedConjugate(twi,conjugator)|
   would make |twi==cs.twistedInvolution(Cartan_nr)| again.
*/

  // Build the parabolic subalgebra:

  cartanclass::InvolutionData id(rd,cs.involutionMatrix(twi));

  PSalgebra result;
  result.cartan = Cartan_nr;

  // Put real simple roots, transformed for original Cartan, into Levi factor
  for (size_t i=0; i <  rd.semisimpleRank(); ++i)
    if (id.real_roots().isMember(rd.simpleRootNbr(i)))
    {
      rootdata::RootNbr alpha=rd.simpleRootNbr(i);
      for (size_t j=conjugator.size(); j-->0; )
	rd.rootReflect(alpha,conjugator[j]);
      result.levi.push_back(alpha);
    }

  // Put pairs of complex positive roots, transformed, into nilpotent radical
  rootdata::RootSet pos_roots=rd.posRootSet();
  for (rootdata::RootSet::iterator i = pos_roots.begin(); i(); ++i)
    if (not id.real_roots().isMember(*i))
    {
      rootdata::RootNbr alpha=*i;
      rootdata::RootNbr beta=id.root_involution()[alpha];
      for (size_t j=conjugator.size(); j-->0; )
      {
	rd.rootReflect(alpha,conjugator[j]);
	rd.rootReflect(beta,conjugator[j]);
      }
      result.nilp.push_back(std::make_pair(alpha,beta));
    }

  return result;

} // theta_stable_parabolic


// Express irreducible K-module as a finite virtual sum of standard ones
CharForm  KHatComputations::character_formula(StandardRepK stdrep) const
{
  const cartanset::CartanClassSet& cs=d_G->cartanClasses();
  weyl::WeylWord conjugator;

  // Get theta stable parabolic subalgebra

  PSalgebra ps = theta_stable_parabolic(conjugator,cs,stdrep.d_cartan);

  // We have the capability to distinguish between compact and non-compact
  // imaginary roots. This information is provided by the fiber group.
  // However, we do not do anything with that infomation yet.

  // Asumme that $u$ is the nilpotent part of |ps| in the Levi deompostion.
  // We will associate to stdrep, $n$ new standardreps, where $n$ is the
  // number of non-empty subsets $A$ of the non-compact of $u$. Each new
  // standardrep will have parameter lambda+sum(of roots in A)

  // For now we will treat the imaginary roots as if they were non compact

  typedef std::pair<atlas::rootdata::RootNbr,atlas::rootdata::RootNbr> rpair;
  std::vector<rpair> rpairlist;
  {
    // Get list of roots in u \cap pc for now we take unique pairs of roots.
    // Later the non compact imaginary roots will be obtained from the fiber

    std::set<rpair> rpairset;

    // make pairs unique by sorting the two elements, then insert into set
    for ( size_t k = 0; k!=ps.nilp.size(); ++k)
    {
      rpair rp = ps.nilp[k];
      if (rp.first > rp.second) std::swap(rp.first,rp.second);
      rpairset.insert(rp);
    }
    rpairlist.assign(rpairset.begin(),rpairset.end());
  }
  size_t u_size = rpairlist.size();

  unsigned long nsubset=1<<u_size; // number of subsets

  normalize(stdrep);

  Char multmap;
  multmap.insert(std::make_pair(stdrep,1));//this handles the empty subset

  latticetypes::LatticeMatrix P=d_data[stdrep.d_cartan].freeProjector;

  for (unsigned long i=1; i<nsubset; ++i) // bits of |i| determine the subset
  {
    bitset::RankFlags subset(i);

    latticetypes::Weight lambda = stdrep.d_lambda.first;
    for (bitset::RankFlags::iterator j =subset.begin(); j(); j++)
      if (rpairlist[*j].first == rpairlist[*j].second) // imaginary root
      { // nothing implemented yet for imaginary roots
      }
      else // complex pair of roots
      {
	latticetypes::Weight mu=
	  P.apply(cs.rootDatum().root(rpairlist[*j].first));
	using latticetypes::operator+=;
	lambda +=mu; // add the restriction of the first root to lambda
      }


    StandardRepK new_rep = stdrep;

    new_rep.d_lambda.first = lambda; // replace |lambda| by modified version
    new_rep.d_status.reset(new_rep.isNormal());

    normalize(new_rep);

    // now add $(-1)^{\#S}$ to the coefficient of |new_rep| in |multmap|
    long sign=subset.count()%2 == 0 ? 1 : -1;
    std::pair<Char::iterator,bool> p=
      multmap.insert(std::make_pair(new_rep,sign));
    if (not p.second) p.first->second += sign;

  } // for (subsets)

  // finally remove items with multiplicity 0

  for (Char::iterator iter = multmap.begin();iter !=multmap.end();)
    if (iter->second == 0)
      multmap.erase(iter++); // must take care to do ++ before erase
    else
      ++iter;

  return std::make_pair(stdrep, multmap);
} // character_formula



/* convert a system of equations into a list, adding equations for all terms
   recursively (they are generated by |character_formula|), up to the given
   bound (everything with |product_sumposroots(...)> bound| is pruned away).
   It is assumed that |character_formula| will ensure all |StandardRepK|s in
   left and right hand sides are normalized, so that there is no risk of
   trying to add a formula for one term but getting one for another.
 */
std::vector<CharForm>
KHatComputations::saturate(std::set<CharForm> system,
			   atlas::latticetypes::LatticeCoeff bound) const
{
  std::set<StandardRepK> lhs; // left hand sides of all formulae seen so far
  for (std::set<CharForm>::iterator it=system.begin(); it!=system.end(); ++it)
    lhs.insert(it->first);

  std::vector<CharForm> result;

  while (not system.empty())
  {
    std::set<CharForm>::iterator current=system.begin(); // choose one
    const CharForm& cf=*current;           // this is an unprocessed formula
    if (product_sumposroots(cf.first) <= bound) // ignore if out of bounds
    {
      result.push_back(CharForm(cf.first,Char())); // start with empty rhs
      Char& rhs=result.back().second;
      for (Char::const_iterator
	     term=cf.second.begin(); term!=cf.second.end(); ++term)
	if (product_sumposroots(term->first) <= bound)
	{
	  rhs.insert(*term);
	  if (lhs.count(term->first)==0) // no formula for this term seen yet
	  {
	    lhs.insert(term->first);
	    system.insert(character_formula(term->first));
	  }
	}

    }
    system.erase(current); // we are done with this formula
  }

  return result;
} // saturate


// **************   manipulators **********************

void KHatComputations::go (const kgb::KGB& kgb)
{
  std::set<StandardRepK> stdRset; // set of standard representations

  rootdata::RootDatum rd = d_G->rootDatum();

  latticetypes::Weight tr=rd.twoRho();
  latticetypes::Weight cotr =rd.dual_twoRho();

  atlas::latticetypes::LatticeCoeff
    bound = atlas::latticetypes::scalarProduct (tr, cotr);
  std::cout << " bound = " << bound << std::endl;

  std::set<CharForm> system;
  for (kgb::KGBElt z =0; z!=kgb.size(); ++z)
  {
    StandardRepK stdrpk=KGB_elt_rep(z);
    //  normalize(*stdrpk);

    CharForm map = character_formula(stdrpk);

    system.insert(map);

  } // for (z)

  std::cout << " size of initial formula set = " << system.size()
	    << std::endl;


  atlas::matrix::Matrix<CharCoeff> m = makeMULTmatrix(system, bound);

} // go



// the following functions have been largely commented out to allow compilation

void KHatComputations::makeHSMap(StandardRepK& stdrep)
{
//   using namespace complexredgp;
//   using namespace latticetypes;
//   using namespace matrix;
//   using namespace rootdata;
//   using namespace cartanset;
//   using namespace descents;

//   // Make a dummy standard rep with

//   StandardRepK undefstdrpk;

//   undefstdrpk.d_cartan = ~0ul;

//   // test if standard

//   if (stdrep.isStandard())
//   {

//     HechtSchmid hs(&undefstdrpk,&undefstdrpk,&undefstdrpk);

//     d_standardHechtSchmid.insert
//       (std::pair<StandardRepK,HechtSchmid>(stdrep,hs));
//   }

//   else
//   {
//     RootDatum rd = d_G->rootDatum();

//     LatticeMatrix theta = d_G->cartan(stdrep.d_cartan).involution();
//     const weyl::TwistedInvolution w =
//       d_G->cartanClasses().twistedInvolution(stdrep.d_cartan);
//     cartanclass::Fiber f = cartanclass::Fiber(rd, theta);
//     atlas::rootdata::RootList rl = f.simpleImaginary();
//     const CartanClassSet& tw = d_G->cartanClasses();

//     Weight lambda = stdrep.d_lambda.first;

//     // make a list of noncompact simple imaginary roots

//     // I am not using isNoncompact yet() because it requires solving systems of
//     // equations

//     rootdata::RootSet ncpi=tw.noncompactRoots(d_realForm);

    // lift lambda to X^*



//   LatticeMatrix basis = d_basis[stdrep.d_cartan];

//   basis.resize(d_basis[stdrep.d_cartan].numRows(), lambda.size());

//   size_t l = d_basis[stdrep.d_cartan].numColumns() - lambda.size();
//   for (size_t k = 0 ; k != basis.numRows(); ++k)
//     for (size_t j  = l; j != d_basis[stdrep.d_cartan].numColumns() ; ++j)
//	basis(k,j-l)=  d_basis[stdrep.d_cartan](k,j);

//   Weight lifted_lambda(basis.numRows());
//   basis.apply(lifted_lambda,lambda);



//     for ( size_t i = 0; i!=rl.size(); ++i)
//     {
// //  I'll have to do this for the first non compact that satisfies the test below
//       if ( ncpi.isMember(rl[i]) &&
// 	   latticetypes::scalarProduct( rd.coroot(rl[i]),lifted_lambda) < 0)

//       {
// 	StandardRepK s, d1, d2;

// 	s= stdrep;
// 	d1 = stdrep;
// 	d2 = undefstdrpk;



// 	// compute s_i(lambda)

// 	LatticeMatrix q;
// 	rd.rootReflection(q,rd.simpleRootNbr(i));
// 	q.apply(lifted_lambda,lifted_lambda);

// 	Weight mu(d_minusQuotient[stdrep.d_cartan].numRows());
//      // minusQuotient(mu, lifted_lambda, stdrep.d_cartan);
// 	s.d_lambda.first = mu;

// 	// Build the standard reps associated to the new Cartan

// 	// Compute Map Lambda->w*lambda'

// 	latticetypes::LatticeMatrix lbdmp = lambdamap(stdrep.d_cartan,i);
// 	Weight wlbdp(lambda.size()-1);

// 	lbdmp.apply(wlbdp, lambda);
// 	d1.d_lambda.first = wlbdp;

// 	// get the new Cartan

// 	weyl::WeylElt w;
// 	size_t j = tw.cayley(stdrep.d_cartan,i,&w);

// 	LatticeMatrix stheta = d_G->cartan(j).involution();
// 	latticetypes::LatticeMatrix qw;

// 	weyl::WeylWord ww;
// 	tw.weylGroup().out(ww,w);
// 	atlas::rootdata::toMatrix(qw, ww,rd);
// 	qw *= d_G->distinguished();
// 	atlas::matrix::conjugate(stheta,qw); // stheta contains theta_prime

// 	// I may need to change this as I do not fully understand ????


// 	d1.d_cartan = j;
// 	// d_G->cartan(j).involution() = stheta;
// 	atlas::descents::DescentStatus ds;
// 	atlas::descents::DescentStatus::Value v =
// 	  ds.operator[](rd.simpleRootNbr(i));

// 	// if i is of type I d2 = d1
// 	if ( v == atlas::descents::DescentStatus::ImaginaryTypeII ) d2 = d1;

// 	HechtSchmid hs(&s,&d1,&d2);
// 	d_standardHechtSchmid.insert
// 	  (std::pair<StandardRepK,HechtSchmid>(stdrep,hs));

// 	break;

//       }

//     }


//  }
  //std::cout << "HS map size " << d_standardHechtSchmid.size() << std::endl;
}

latticetypes::LatticeMatrix
KHatComputations::lambdamap(size_t cartan_num,size_t i)
{
  using namespace complexredgp;
  using namespace latticetypes;
  using namespace matrix;
  using namespace rootdata;
  using namespace cartanset;



  RootDatum rd = d_G->rootDatum();
  const CartanClassSet& tw = d_G->cartanClasses();
  weyl::WeylElt w;

  latticetypes::LatticeMatrix qw,m;
  weyl::WeylWord ww;
  tw.weylGroup().out(ww,w);
  atlas::rootdata::toMatrix(qw, ww,rd);
  qw*= d_G->cartanClasses().distinguished();

//size_t j = tw.cayley(cartan_num,i,&w);
//qw*=d_basisInverse[cartan_num];
//m = d_basis[j];
//m *=qw; // this matrix is the product of 3 n by n matrices

//size_t r = d_minusQuotient[cartan_num].numRows();
//size_t n = m.numRows();

  latticetypes::LatticeMatrix ldm; //(r-1,r);


//   for ( size_t k = 0; k!=r-1; ++k)
//     for ( size_t l = 0 ; l!=r; ++l)
//       ldm(k,l) = m(n-r+1+k,n-r+l);

  return ldm;

}


// ****************** functions ************************

atlas::matrix::Matrix<CharCoeff>
triangularize (const std::vector<CharForm>& system,
	       std::vector<StandardRepK>& new_order)
{
  std::vector<CharForm> equation(system.begin(),system.end()); // numbering
  size_t n=equation.size();

  atlas::matrix::Matrix<CharCoeff> M(n,n,0);
  graph::OrientedGraph usage(n);
  for (size_t i=0; i<n; ++i) // loop over equations
  {
    size_t n_terms=0;
    for (size_t j=0; j<n; ++j) // loop over left hand sides
    {
      Char::const_iterator p= equation[i].second.find(equation[j].first);
      if (p!=equation[i].second.end())
      { // |OrientedGraph::cells| put sinks in front, so record edge $j\to i$.
	usage.edgeList(j).push_back(i); M(i,j)=p->second;
	++n_terms;
      }
    }
    if (equation[i].second.size()!=n_terms)
      throw std::runtime_error ("triangularize: system not saturated");
  }

  partition::Partition order; usage.cells(order,NULL);

  new_order.resize(n);
  for (size_t i=0; i<n; ++i)
  {
    if (order.classSize(i)>1)
      throw std::runtime_error ("triangularize: system has cycles");
    new_order[order(i)]=equation[i].first;
  }

  atlas::matrix::Matrix<CharCoeff> result(n,n,0);
  for (size_t j=0; j<n; ++j)
    for (graph::EdgeList::const_iterator
	   p=usage.edgeList(j).begin(); p != usage.edgeList(j).end(); ++p)
      result(order(*p),order(j))=M(*p,j);
  return result;
} // triangularize

matrix::Matrix<CharCoeff> inverse_upper_triangular
  (const atlas::matrix::Matrix<CharCoeff>& U)
{
  size_t n=U.numColumns();
  if (U.numRows()!=n)
    throw std::runtime_error ("invert triangular: matrix is not square");

  matrix::Matrix<CharCoeff> result(n,n,0);

  for (size_t j=0; j<n; ++j)
  {
    if (U(j,j)!=1)
      throw std::runtime_error ("invert triangular: not unitriangular");
    result(j,j)=1;

    for (size_t i=j; i-->0; )
    {
      CharCoeff sum=0;
      for (size_t k=j; k>i; --k)
	sum += U(i,k)*result(k,j);
      result(i,j) = -sum;
    }
  }
  return result;
}

atlas::matrix::Matrix<CharCoeff>
KHatComputations::makeMULTmatrix
 (std::set<CharForm>& column,
  const atlas::latticetypes::LatticeCoeff bound)
{
  using latticetypes::operator+=;

  rootdata::RootDatum rd = d_G->rootDatum();
  latticetypes::Weight tr=rd.twoRho();
  latticetypes::Weight cotr = rd.dual_twoRho();

  std::set <atlas::latticetypes::Weight> lookup;

  //it is assumed that the characters in column.first are different

  for (std::set<CharForm>::iterator i = column.begin(); i != column.end();)
    if (product_sumposroots(i->first) > bound)
      column.erase(i++); // increment iterator before actual erase
    else
    {
      lookup.insert(i->first.d_lambda.first);
      ++i;
    }

  assert (column.size()==lookup.size());


  std::vector<CharForm> system=saturate(column,bound);

  std::vector<StandardRepK> new_order;
  atlas::matrix::Matrix<CharCoeff>  m=triangularize(system,new_order);
  for (std::vector<StandardRepK>::const_iterator
	 it=new_order.begin(); it!=new_order.end(); ++it)
  {
    const latticetypes::Weight& l=it->d_lambda.first;
    basic_io::seqPrint(std::cout,l.begin(),l.end(), ", ", "[", "]\n");
  }

  prettyprint::printMatrix(std::cout<<"Triangular system:\n",m,3);

  matrix::Matrix<CharCoeff>m_inv=inverse_upper_triangular(m);

  prettyprint::printMatrix(std::cout<<"Inverse matrix:\n",m_inv,3);

  return m_inv;

}




} // namespace standardrepk
} // namespace atlas
