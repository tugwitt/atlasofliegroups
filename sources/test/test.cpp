/*
  This is test.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "io.h"    //needed for help commands

#include "test.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

#include "atlas_types.h" // here to preempt double inclusion of _fwd files

#include "free_abelian.h"
#include "permutations.h"
#include "matreduc.h"

#include "dynkin.h"
#include "prerootdata.h"
#include "rootdata.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "realredgp.h"

#include "kgb.h"
#include "blocks.h"
#include "klsupport.h"
#include "kl.h"
#include "standardrepk.h"
#include "repr.h"
#include "ext_block.h"

#include "ioutils.h"
#include "basic_io.h"
#include "hkl.h"
#include "prettyprint.h"
#include "interactive.h"
#include "realform_io.h"
#include "kgb_io.h"
#include "block_io.h"

#include "commands.h"
#include "helpmode.h"
#include "emptymode.h"
#include "mainmode.h"
#include "realmode.h"
#include "blockmode.h"
#include "reprmode.h"

#include "testrun.h"
#include "kltest.h"


namespace atlas {

namespace test {

/*****************************************************************************

  This module contains some commands for testing the program.

******************************************************************************/

namespace {
  // functions for the test commands

  void test_f();

  void roots_rootbasis_f();
  void coroots_rootbasis_f();
  void posroots_rootbasis_f();
  void poscoroots_rootbasis_f();
  void Ktypeform_f();
  void qKtypeform_f();
  void Ktypemat_f();
  void qKtypemat_f();
  void mod_lattice_f();
  void branch_f();
  void qbranch_f();
  void srtest_f();

  void X_f();
  void hblock_f();
  void hkllist_f();
  void hklbasis_f();
  void hkltest_f();

/*
  For convenience, the "test" command is added to the mode that is flagged by
  the testMode constant defined here; therefore "test" appears (conditionally)
  in every template instance of |addTestCommands| and of |addTestHelp| below.
  Set this constant according to the requirements of the |test_f| function.
*/
  enum TestMode {EmptyMode, MainMode, RealMode, BlockMode, ReprMode,
		 numTestMode};
  const TestMode testMode = ReprMode; // currently does test of extended block

  // utilities
  const RootDatum& currentRootDatum();

} // |namespace|

/*****************************************************************************

        Chapter I -- Functions declared by test.h

  This section defines the functions declared in test.h :

    - addTestCommands() : adds the test commands to the main command tree;

******************************************************************************/

  const char* test_tag = "test command (for development only)";


/* |addTestCommands| is a template only so that its declaration is shorter;
   its defining instances are defined just as overloaded functions would,
   since each of them needs to test a specific value of |testMode|.
*/


// Add to the empty mode node the test commands that require that mode.
template<>
void addTestCommands<commands::EmptymodeTag> (commands::CommandNode& mode)
{
  if (testMode == EmptyMode)
    mode.add("test",test_f,test_tag);
}


// Add to the main mode node the test commands that require that mode.
template<>
void addTestCommands<commands::MainmodeTag> (commands::CommandNode& mode)
{
  mode.add("roots_rootbasis",roots_rootbasis_f,
	   "outputs the roots in the simple root basis",commands::std_help);
  mode.add("posroots_rootbasis",posroots_rootbasis_f,
	   "outputs the positive roots in the simple root basis",
	   commands::std_help);
  mode.add("coroots_rootbasis",coroots_rootbasis_f,
	   "outputs the coroots in the simple coroot basis",
	   commands::std_help);
  mode.add("poscoroots_rootbasis",poscoroots_rootbasis_f,
	   "outputs the positive coroots in the simple coroot basis",
	   commands::std_help);

  mode.add("X",X_f,"prints union of K\\G/B for real forms in inner class",
	   commands::std_help);

  if (testMode == MainMode)
    mode.add("test",test_f,test_tag);

}

// Add to the real mode node the test commands that require that mode.
template<>
void addTestCommands<commands::RealmodeTag> (commands::CommandNode& mode)
{
  mode.add("Ktypeform",Ktypeform_f,
	   "computes formula for a K-type",commands::std_help);
  mode.add("qKtypeform",qKtypeform_f,
	   "q version of Ktypeform",commands::use_tag);
  mode.add("Ktypemat",Ktypemat_f,
	   "computes matrix relating K-types and standard modules",
	   commands::std_help);
  mode.add("qKtypemat",qKtypemat_f,"q version of Ktypemat",commands::use_tag);
  mode.add("mod_lattice",mod_lattice_f,
	   "gives a basis of quotient of character lattice",commands::std_help);
  mode.add("branch",branch_f,
	   "computes restriction of representation to K",commands::std_help);
  mode.add("qbranch",qbranch_f,"q version of branch");
  mode.add("srtest",srtest_f,
	   "gives information about a representation",commands::std_help);

  if (testMode == RealMode)
    mode.add("test",test_f,test_tag);

}


// Add to the block mode the test commands that require that mode.
template<>
void addTestCommands<commands::BlockmodeTag> (commands::CommandNode& mode)
{
  mode.add("hblock",hblock_f,
	   "hermitian variant of block command for complex groups",
	   commands::use_tag);
  mode.add("hkllist",hkllist_f,
	   "hermitian variant of kllist command for complex groups",
	   commands::use_tag);
  mode.add("hklbasis",hklbasis_f,
	   "hermitian variant of klbasis command for complex groups",
	   commands::use_tag);
  mode.add("hkltest",hkltest_f,
	   "hermitian variant of kltest command for complex groups",
	   commands::use_tag);

  if (testMode == BlockMode)
    mode.add("test",test_f,test_tag);
}

// Add to the repr mode the test commands that require that mode.
template<>
void addTestCommands<commands::ReprmodeTag> (commands::CommandNode& mode)
{
  if (testMode == ReprMode)
    mode.add("test",test_f,test_tag);

  // add additional commands here :
}


/*****************************************************************************

        Chapter II -- Functions for the test commands

******************************************************************************/

namespace {


// Empty mode functions


// Main mode functions

// Print the roots in the simple root coordinates.
void roots_rootbasis_f()
{
  const RootSystem& rs =  commands::currentComplexGroup().rootSystem();
  ioutils::OutputFile file;

  for (RootNbr i=0; i<rs.numRoots(); ++i)
    prettyprint::printInRootBasis(file,i,rs) << std::endl;
}

// Print the positive roots in the simple root coordinates.
void posroots_rootbasis_f()

{
  const RootSystem& rs = commands::currentComplexGroup().rootSystem();

  ioutils::OutputFile file;
  prettyprint::printInRootBasis(file,rs.posRootSet(),rs);
}

// Print the coroots in the simple coroot coordinates.
void coroots_rootbasis_f()
{
  const RootSystem rs (commands::currentComplexGroup().dualRootSystem());

  ioutils::OutputFile file;
  for (RootNbr i=0; i<rs.numRoots(); ++i)
    prettyprint::printInRootBasis(file,i,rs) << std::endl;

}

// Print the positive coroots in the simple coroot coordinates.
void poscoroots_rootbasis_f()
{
  const RootSystem rs (commands::currentComplexGroup().dualRootSystem());

  ioutils::OutputFile file;
  prettyprint::printInRootBasis(file,rs.posRootSet(),rs);
}


void X_f()
{
  ComplexReductiveGroup& G=commands::currentComplexGroup();
  kgb::global_KGB kgb(G); // build global Tits group, "all" square classes
  ioutils::OutputFile f;
  kgb_io::print_X(f,kgb);
}


// Real mode functions

void Ktypeform_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::KhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  ioutils::OutputFile f;

  standardrepk::CharForm kf= khc.K_type_formula(sr);

  khc.print(f << "K-type formula for mu(",kf.first) << "):\n";
  {
    std::ostringstream s; khc.print(s,kf.second);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
  }

  for (size_t i=0; i<khc.nr_reps(); ++i)
    khc.print(f << 'R' << i << ": ",khc.rep_no(i))
      << ", height: " << khc.height(i) << std::endl;


  standardrepk::combination sum(khc.height_order());

  for (standardrepk::Char::const_iterator
	 it=kf.second.begin(); it!=kf.second.end(); ++it)
  {
#ifdef VERBOSE
    khc.print(f,it->first) << " has height " << khc.height(it->first)
				   << std::endl;
    size_t old_size=khc.nr_reps();
#endif
    standardrepk::combination st=khc.standardize(it->first);
#ifdef VERBOSE
    for (size_t i=old_size; i<khc.nr_reps(); ++i)
      khc.print(f << 'R' << i << ": ",khc.rep_no(i))
        << ", height: " << khc.height(i) << std::endl;

    std::ostringstream s; khc.print(s,it->first) << " = ";
    khc.print(s,st,true);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
#endif
    sum.add_multiple(st,it->second);
  }

  f << "Converted to Standard normal final limit form:\n";
  {
    std::ostringstream s; khc.print(s,sum);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
  }
} // |Kypeform_f|

void qKtypeform_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::qKhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  ioutils::OutputFile f;

  standardrepk::q_CharForm kf= khc.q_K_type_formula(sr);

  khc.print(f << "q-K-type formula for mu(",kf.first) << "):\n";
  {
    std::ostringstream s; khc.print(s,kf.second);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
  }

  standardrepk::q_combin sum(khc.height_order());

  for (standardrepk::q_Char::const_iterator
	 it=kf.second.begin(); it!=kf.second.end(); ++it)
  {
#ifdef VERBOSE
    khc.print(f,it->first) << " has height " << khc.height(it->first)
				   << std::endl;
    size_t old_size=khc.nr_reps();
#endif
    standardrepk::q_combin st=khc.standardize(it->first);
#ifdef VERBOSE
    for (size_t i=old_size; i<khc.nr_reps(); ++i)
      khc.print(f << 'R' << i << ": ",khc.rep_no(i))
        << ", height: " << khc.height(i) << std::endl;

    std::ostringstream s; khc.print(s,it->first) << " = ";
    khc.print(s,st,true);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
#endif
    sum.add_multiple(st,it->second);
  }

  f << "Converted to Standard normal final limit form:\n";
  {
    std::ostringstream s; khc.print(s,sum);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
  }
} // |qKypeform_f|

void Ktypemat_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::KhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);
  khc.normalize(sr);

  {
    size_t witness;
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  standardrepk::combination c=khc.standardize(sr);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",sr) << " is zero.\n";
    return;
  }

  assert(c.size()==1);
  assert(khc.rep_no(c.begin()->first)==sr);

  khc.print(std::cout << "Height of representation ",sr) << " is "
    << khc.height(c.begin()->first) << ".\n";
  unsigned long bound=
    interactive::get_bounded_int(interactive::common_input(),
				 "Give height bound: ",
				 9999);

  ioutils::OutputFile f;

  std::set<standardrepk::equation> singleton;
  singleton.insert(khc.mu_equation(c.begin()->first,bound));

  {
    standardrepk::equation init=*singleton.begin();
    khc.print(f << "Initial formula: mu(",khc.rep_no(init.first))
      << ") =\n";
    {
      std::ostringstream s; khc.print(s,init.second);
      ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
    }
  }
  std::vector<standardrepk::seq_no> new_order;

#ifdef VERBOSE
  matrix::Matrix_base<standardrepk::CharCoeff> m;
  matrix::Matrix_base<standardrepk::CharCoeff> ktypemat =
    khc.K_type_matrix(singleton,bound,new_order,&m);
#else
  matrix::Matrix_base<standardrepk::CharCoeff> ktypemat =
    khc.K_type_matrix(singleton,bound,new_order,NULL);
#endif

  std::cout << "Ordering of representations/K-types:\n";
  for (std::vector<standardrepk::seq_no>::const_iterator
	 it=new_order.begin(); it!=new_order.end(); ++it)
    khc.print(std::cout,khc.rep_no(*it)) << ", height " << khc.height(*it)
       << std::endl;

#ifdef VERBOSE
  prettyprint::printMatrix(std::cout<<"Triangular system:\n",m,3);
#endif

  prettyprint::printMatrix(f<<"Matrix of K-type multiplicites:\n",ktypemat,3);

} // |Ktypemat_f|

void qKtypemat_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::qKhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isNormal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not normal, as witnessed by coroot sum "
	<< khc.info(sr.Cartan()).coroot_sum(witness)
	<< ".\n";
      return;
    }
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  standardrepk::q_combin c=khc.standardize(sr);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",sr) << " is zero.\n";
    return;
  }

  assert(c.size()==1);
  assert(khc.rep_no(c.begin()->first)==sr);

  khc.print(std::cout << "Height of representation ",sr) << " is "
		      << khc.height(c.begin()->first) << ".\n";
  unsigned long bound=
    interactive::get_bounded_int(interactive::common_input(),
				 "Give height bound: ",
				 9999);

  ioutils::OutputFile f;

  std::set<standardrepk::q_equation> singleton;
  singleton.insert(khc.mu_equation(c.begin()->first,bound));

  {
    standardrepk::q_equation init=*singleton.begin();
    khc.print(f << "Initial formula: mu(",khc.rep_no(init.first))
      << ") =\n";
    {
      std::ostringstream s; khc.print(s,init.second);
      ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
    }
  }

#ifdef VERBOSE

  std::vector<standardrepk::q_equation> system =
    khc.saturate(singleton,bound);
  std::cout << "System of equations:\n";
  for (size_t i=0; i<system.size(); ++i)
  {
    const standardrepk::q_equation& si=system[i];
    khc.print(std::cout<< si.first << ' ',khc.rep_no(si.first))
      << " [" << khc.height(si.first) << "]\n     ";
    for (standardrepk::q_combin::const_iterator
	   it=si.second.begin(); it!=si.second.end(); ++it)
      std::cout << '+' << it->second << "*I(" << it->first << ')';
    std::cout << std::endl;
  }
#endif

  std::vector<standardrepk::seq_no> new_order;

#ifdef VERBOSE
  matrix::Matrix_base<standardrepk::q_CharCoeff> m;
  matrix::Matrix_base<standardrepk::q_CharCoeff> ktypemat =
    khc.K_type_matrix(singleton,bound,new_order,&m);
#else
  matrix::Matrix_base<standardrepk::q_CharCoeff> ktypemat =
    khc.K_type_matrix(singleton,bound,new_order,NULL);
#endif

  f << "Ordering of representations/K-types:\n";
  for (std::vector<standardrepk::seq_no>::const_iterator
	 it=new_order.begin(); it!=new_order.end(); ++it)
    khc.print(f,khc.rep_no(*it)) << ", height " << khc.height(*it)
				 << std::endl;

#ifdef VERBOSE
  prettyprint::printMatrix(std::cout<<"Triangular system:\n",m,3);
#endif

  prettyprint::printMatrix(f<<"Matrix of K-type multiplicites:\n",ktypemat,3);
} // |qKtypemat_f|

void mod_lattice_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();

  unsigned long cn=interactive::get_Cartan_class(G.Cartan_set());

  WeightInvolution q = G.cartan(cn).involution();
  for (size_t j = 0; j<q.numRows(); ++j)
    q(j,j) -= 1;

  CoeffList factor;
  int_Matrix b = matreduc::adapted_basis(q,factor);

  RankFlags units, doubles;
  unsigned n1=0,n2=0;

  for (size_t i=0; i<factor.size(); ++i)
    if (factor[i]==1)
      units.set(i),++n1;
    else if (factor[i]==2)
      doubles.set(i),++n2;

   std::cout << "At Cartan class " << cn;
   if (n1+n2==0)
     std::cout << " weights are used unchanged";
   else
   {
     std::cout << " weights are modulo";
     if (n1>0)
     {
       std::cout << " multiples of ";
       for (RankFlags::iterator it=units.begin(); it(); ++it,--n1)
	 std::cout << b.column(*it) << (n1>2 ? ", " : n1>1 ? ", and " : "");
       if (n2>0)
	 std::cout << " and";
     }
     if (n2>0)
     {
       std::cout << " even multiples of ";
       for (RankFlags::iterator it=doubles.begin(); it(); ++it,--n2)
	 std::cout << b.column(*it) << (n2>2 ? ", " : n2>1 ? ", and " : "");
     }
   }
   std::cout << ".\n";

} // |mod_lattice_f|

void branch_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::KhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  standardrepk::combination c=khc.standardize(sr);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",sr) << " is zero.\n";
    return;
  }

  assert(c.size()==1 and khc.rep_no(c.begin()->first)==sr);

  khc.print(std::cout << "Height of representation ",sr) << " is "
    << khc.height(c.begin()->first) << ".\n";
  unsigned long bound=
    interactive::get_bounded_int(interactive::common_input(),
				 "Give height bound: ",
				 9999);

  ioutils::OutputFile f;

  standardrepk::combination result=khc.branch(c.begin()->first,bound);

  {
    std::ostringstream s; khc.print(s,result);
    ioutils::foldLine(f,s.str(),"+","",1) << std::endl;
  }

} // |branch_f|

void qbranch_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::qKhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  standardrepk::q_combin c=khc.standardize(sr);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",sr) << " is zero.\n";
    return;
  }

  assert(c.size()==1 and khc.rep_no(c.begin()->first)==sr);

  khc.print(std::cout << "Height of representation ",sr) << " is "
    << khc.height(c.begin()->first) << ".\n";
  unsigned long bound=
    interactive::get_bounded_int(interactive::common_input(),
				 "Give height bound: ",
				 9999);

  ioutils::OutputFile f;

  standardrepk::q_combin result=khc.branch(c.begin()->first,bound);

  {
    std::ostringstream s; khc.print(s,result);
    ioutils::foldLine(f,s.str(),"+","",1) << std::endl;
  }

} // |qbranch_f|


/*
  Function invoked by the "srtest" command.
*/
void srtest_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();
  const KGB& kgb = G.kgb();

  unsigned long x=interactive::get_bounded_int
    (interactive::sr_input(),"Choose KGB element: ",G.kgb().size());

  prettyprint::printVector(std::cout<<"2rho = ",G.rootDatum().twoRho())
    << std::endl;

  Weight lambda=
    interactive::get_weight(interactive::sr_input(),
			    "Give lambda-rho: ",
			    G.rank());
  standardrepk::KhatContext khc(G);

  StandardRepK sr=khc.std_rep_rho_plus(lambda,kgb.titsElt(x));

  (lambda *= 2) += G.rootDatum().twoRho();
  prettyprint::printVector(std::cout << "Weight (1/2)",lambda);
  prettyprint::printVector(std::cout << " converted to (1/2)",khc.lift(sr));

  const TwistedInvolution& canonical =
    G.complexGroup().twistedInvolution(sr.Cartan());
  if (kgb.involution(x)!=canonical)
    prettyprint::printWeylElt(std::cout << " at involution ",
			      canonical, G.weylGroup());
  std::cout << "\nHeight is " << khc.height(sr) << std::endl;

  khc.go(sr);
}

bool examine(RealReductiveGroup& G)
{
  const WeylGroup& W = G.weylGroup();
  const KGB& kgb=G.kgb();
  size_t l = W.length(kgb.involution(0)),t;
  for (size_t i=1; i<kgb.size(); ++i)
    if ((t=W.length(kgb.involution(i)))<l)
      return false;
    else
      l=t;
  return true;
}



TorusElement torus_part
  (const RootDatum& rd,
   const WeightInvolution& theta,
   const RatWeight& lambda, // discrete parameter
   const RatWeight& gamma // infinitesimal char
  )
{
  InvolutionData id(rd,theta);
  Weight cumul(rd.rank(),0);
  arithmetic::Numer_t n=gamma.denominator();
  const Ratvec_Numer_t& v=gamma.numerator();
  const RootNbrSet pos_real = id.real_roots() & rd.posRootSet();
  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
    if (rd.coroot(*it).dot(v) %n !=0) // nonintegral
      cumul+=rd.root(*it);
  // now |cumul| is $2\rho_\Re(G)-2\rho_\Re(G(\gamma))$

  return y_values::exp_pi(gamma-lambda+RatWeight(cumul,2));
}

void test_f()
{
} // |test_f|


// Block mode functions

void hblock_f()
{
  try
  {
    blocks::hBlock block = hBlock(commands::currentBlock());

    ioutils::OutputFile f;
    block_io::printHBlockD(f,block);

    Block* block_pointer = new Block
      (Block::build(commands::currentRealGroup(),
		    commands::currentDualRealGroup()));
    kl::KLContext* klc_pointer=new kl::KLContext(*block_pointer);
    klc_pointer->fill();

    size_t hbsize = block.hsize();
    size_t klsize = klc_pointer->size();
    std::vector<bool> flags(klsize,false);
    for (size_t i=0; i<hbsize; i++) {
      BlockElt j = block.hfixed(i);
      flags[j] = true;
    }

    size_t count = 0;
    int width = ioutils::digits(klc_pointer->size()-1,10ul);
    int tab = 2;

    for (size_t y = 0; y < klc_pointer->size(); ++y)
      if (flags[y])
      {
	f << std::setw(width) << y << ": ";
	bool first = true;

	for (size_t x = 0; x <= y; ++x)
	  if (flags[x])
	  {
	    const kl::KLPol& pol = klc_pointer->klPol(x,y);
	    if (pol.isZero())
	      continue;
	    if (first)
	    {
	      f << std::setw(width) << x << ": ";
	      first = false;
	    }
	    else
	    {
	      f << std::setw(width+tab)<< ""
		<< std::setw(width) << x << ": ";
	    }
	    pol.print(f,"q");
	    f << "" << std::endl;
	    ++count;
	  }

	f << "" << std::endl;
      }

    /*
      kl::KLContext klc(block);
      klc.fill(z,false);

      typedef Polynomial<int> Poly;
      typedef std::map<BlockElt,Poly> map_type;
      map_type acc;
      unsigned int parity = block.length(z)%2;
      for (size_t x = 0; x <= z; ++x)
      {
      const kl::KLPol& pol = klc.klPol(x,z);
      if (not pol.isZero())
      {
      Poly p(pol); // convert
      if (block.length(x)%2!=parity)
      p*=-1;
      BlockEltList nb=block.nonzeros_below(x);
      for (size_t i=0; i<nb.size(); ++i)
      {
      std::pair<map_type::iterator,bool> trial =
      acc.insert(std::make_pair(nb[i],p));
      if (not trial.second) // failed to create a new entry
      trial.first->second += p;
      } // |for (i)| in |nb|
      } // |if(pol!=0)|
      } // |for (x<=z)|


      f << (block.singular_simple_roots().any() ? "(cumulated) " : "")
      << "KL polynomials (-1)^{l(" << z << ")-l(x)}*P_{x," << z << "}:\n";
      int width = ioutils::digits(z,10ul);
      for (map_type::const_iterator it=acc.begin(); it!=acc.end(); ++it)
      {
      BlockElt x = it->first;
      const Poly& pol = it->second;
      if (not pol.isZero())
      {
      f << std::setw(width) << x << ": ";
      pol.print(f,"q") << std::endl;
      }
      }
    */
  }
  catch (error::MemoryOverflow& e)
  {
    e("error: memory overflow");
  }
  catch (error::InputError& e)
  {
    e("aborted");
  }
  catch (std::exception& e)
  {
    std::cerr << "error occurrend: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << std::endl << "unidentified error occurred" << std::endl;

  }

} // |hblock_f|

void hkllist_f()
{
  ioutils::OutputFile file;
  std::ostream& strm = file;

  blocks::hBlock block = hBlock(commands::currentBlock());

  std::cerr << "Computing twisted KL polynomials ... ";
  kl::hKLContext* hklc = new kl::hKLContext(block);
  hklc->fill();
  std::cerr << "Done." << std::endl;

  // sort the polynomials before printing
  const kl::hKLStore& plist = hklc->getPolyList();
  std::vector<kl::hKLPol> splist;
  size_t psize = plist.size();
  for (kl::KLIndex i=1; i<psize; ++i)
    splist.push_back(plist[i]);
  strm << "Printing " << psize-1 << " twisted KL polynomials:" << std::endl;

  // sort
  std::sort(splist.begin(),splist.end(),polynomials::compare<kl::hKLCoeff>);

  // write the polynomials
  for (size_t i=0; i<psize-1; ++i)
    splist[i].print(strm,"q") << std::endl;
}

void hklbasis_f()
{
  ioutils::OutputFile file;
  std::ostream& strm = file;

  blocks::hBlock block = hBlock(commands::currentBlock());

  std::cerr << "Computing twisted KL polynomials ... ";
  kl::hKLContext* hklc = new kl::hKLContext(block);
  hklc->fill();
  std::cerr << "Done." << std::endl;

  // print the full basis
  size_t hsize = hklc->getSize();
  int width = ioutils::digits(hsize-1,10ul);
  strm << "Printing " << hsize << " elements with "
       << hklc->getPolyList().size()-1 << " distinct polynomials:"
       << std::endl;
  for (size_t y=0; y<hsize; ++y)
  {
    //for (size_t y=hsize; y-->0;) {
    strm << std::setw(width) << y << ": ";
    bool first = true;

    for (size_t x = 0; x <= y; ++x)
    {
      const kl::hKLPol& pol = hklc->getPoly(x,y);
      if (pol.isZero()) continue;
      if (first)
      {
	strm << std::setw(width) << x << ": ";
	first = false;
      }
      else
	strm << std::setw(width+2)<< "" << std::setw(width) << x << ": ";
      pol.print(strm,"q");
      strm << std::endl;
    }

    strm << std::endl;
  }
}

void hkltest_f()
{
  ioutils::OutputFile file;
  std::ostream& strm = file;

  blocks::hBlock hblock(commands::currentBlock());

  std::cerr << "Computing KL polynomials ... ";
  kl::KLContext klc(hblock.base());
  klc.fill();
  std::cerr << "Done." << std::endl;

  std::cerr << "Computing twisted KL polynomials ... ";
  kl::hKLContext hklc(hblock);
  hklc.fill();
  std::cerr << "Done." << std::endl;

  std::cerr << "Checking polynomials ... ";
  // compare the two sets of polynomials and test for obvious errors
  size_t count = 0;
  size_t hbsize = hblock.hsize();
  for (BlockElt y=0; y<hbsize; ++y)
  {
    BlockElt by = hblock.hfixed(y);
    for (BlockElt x=0; x<=y; ++x)
    {
      BlockElt bx = hblock.hfixed(x);
      const kl::KLPol& pol = klc.klPol(bx,by);
      const kl::hKLPol& hpol = hklc.getPoly(x,y);
      if (pol != hpol) {
	size_t d1 = pol.degree();
	size_t d2 = hpol.degree();
	size_t top = (d1 >= d2) ? d1 : d2;

	// check the coefficients
	for (size_t i=0; i<=top; ++i)
	{
	  kl::hKLCoeff c1 = (i<=d1) ? pol[i] : 0;
	  kl::hKLCoeff c2 = (i<=d2) ? hpol[i] : 0;
	  if (c2 < 0)
	    c2 = -c2;
	  assert (c2 <= c1);
	  assert ((c1-c2) % 2 ==0);
	}
      }
    }
  }
  std::cerr << "Passed!" << std::endl;

  int tab = 2;
  int width = ioutils::digits(hbsize-1,10ul);
  strm << "Printing " << hbsize << " elements with "
       << hklc.getPolyList().size()-1 << " distinct polynomials:"
       << std::endl;
  for (BlockElt y=0; y<hbsize; y++)
  {
    bool first = true;
    BlockElt by = hblock.hfixed(y);
    strm << std::setw(width) << y << ": ";

    for (BlockElt x=0; x<=y; x++)
    {
      BlockElt bx = hblock.hfixed(x);
      const kl::KLPol& pol = klc.klPol(bx,by);
      const kl::hKLPol& hpol = hklc.getPoly(x,y);

      if (pol.isZero() and hpol.isZero())
	continue;
      if (first)
      {
	strm << std::setw(width) << x << ": ";
	pol.print(strm,"q");
	if (hpol != pol)
	{
	  strm << std::endl;
	  strm << std::setw(2*width+2*tab) << "";
	  hpol.print(strm,"q");
	}
	first = false;
      }
      else
      {
	strm << std::setw(width+tab) << "" << std::setw(width) << x << ": ";
	pol.print(strm,"q");
	if (hpol != pol)
	{
	  strm << std::endl;
	  strm << std::setw(2*width+2*tab) << "";
	  hpol.print(strm,"q");
	}
      }
      strm << "" << std::endl;
      ++count;
    }

    strm << "" << std::endl;
  }

} // |hkltest_f|



} // |namespace|

} // |namespace test|

} // |namespace atlas|

