/*
  This is test.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "test.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include "basic_io.h"
#include "bitmap.h"
#include "block_io.h"
#include "blocks.h"
#include "bruhat.h"
#include "cartan_io.h"
#include "cartanset.h"
#include "commands.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "dynkin.h"
#include "error.h"
#include "gradings.h"
#include "input.h"
#include "interactive.h"
#include "involutions.h"
#include "io.h"
#include "ioutils.h"
#include "kl.h"
#include "kl_io.h"
#include "klsupport.h"
#include "kltest.h"
#include "kgb.h"
#include "kgb_io.h"
#include "lattice.h"
#include "latticetypes.h"
#include "poset_io.h"
#include "prettyprint.h"
#include "realform.h"
#include "realform_io.h"
#include "realredgp_io.h"
#include "realweyl.h"
#include "realweyl_io.h"
#include "rootdata.h"
#include "size.h"
#include "smithnormal.h"
#include "tits.h"
#include "tori.h"
#include "weyl.h"
#include "wgraph.h"
#include "wgraph_io.h"

#include "helpmode.h"
#include "emptymode.h"
#include "mainmode.h"
#include "realmode.h"

#include "testprint.h"
#include "testrun.h"
#include "filekl.h"

/*****************************************************************************

  This module contains some commands for testing the program.

******************************************************************************/

namespace atlas {

namespace {
  using namespace test;

  // functions for the test commands

  void cmatrix_f();
  void corder_f();
  void components_f();
  void coroots_rootbasis_f();
  void poscoroots_rootbasis_f();
  void posroots_rootbasis_f();
  void roots_rootbasis_f();
  void rootdatum_f();
  void kgborder_f();
  void primkl_f();
template<bool small>
  void kgb_f();
template<bool small>
  void block_f();
template<bool small>
  void dual_kgb_f();
template<bool small>
  void dual_block_f();
  void dual_map_f();
  void blockd_f();
  void blockorder_f();
  void blocku_f();
  void blockstabilizer_f();
  void klbasis_f();
  void kllist_f();
  void wcells_f();
  void wgraph_f();
  void klwrite_f();
  void blockwrite_f();
  void extract_graph_f();
  void extract_cells_f();
  void test_f();

  // help functions

  void block_h();
  void blockd_h();
  void blockorder_h();
  void blocku_h();
  void cmatrix_h();
  void kgborder_h();
  void primkl_h();
  void kgb_h();
  void klbasis_h();
  void kllist_h();
  void klwrite_h();
  void blockwrite_h();
  void wcells_h();
  void wgraph_h();

  // tags

  const char* test_tag = "(test command)";
  const char* block_tag = "prints all the representations in a block";
  const char* blocku_tag =
   "prints the unitary representations in the block at rho";
  const char* cmatrix_tag = "prints the Cartan matrix";
  const char* kgb_tag = "prints the orbits of K on G/B";
  const char* dual_kgb_tag = "prints the KGB data for a dual real form";
  const char* dual_block_tag = "prints a block for the dual group";
  const char* dual_map_tag = "prints a map from block to its dual block";
  const char* klbasis_tag = "prints the KL basis for the Hecke module";
  const char* kllist_tag = "prints the list of distinct KL polynomials";
  const char* klwrite_tag = "writes the KL polynomials to disk";
  const char* blockwrite_tag = "writes the block information to disk";
  const char* wcells_tag = "prints the Kazhdan-Lusztig cells for the block";
  const char* wgraph_tag = "prints the W-graph for the block";
  const char* extract_graph_tag =
   "reads block and KL binary files and prints W-graph";
  const char* extract_cells_tag =
   "reads block and KL binary files and prints W-cells";

/*
  For convenience, the "test" command is added to the mode that is flagged by
  the testMode constant defined here; therefore "test" appears (conditionally)
  in every template instance of |addTestCommands| and of |addTestHelp| below.
  Set this constant according to the requirements of the |test_f| function.
*/
  enum TestMode {EmptyMode, MainMode, RealMode, numTestMode};
  const TestMode testMode = RealMode; // currently does a KL matrix test

  // utilities
  const rootdata::RootDatum& currentRootDatum();

}

/*****************************************************************************

        Chapter I -- Functions declared by test.h

  This section defines the functions declared in test.h :

    - addTestCommands() : adds the test commands to the main command
      tree;
    - addTestHelp() : adds help functionality;

******************************************************************************/

namespace test {

/* |addTestCommands| is a template only so that its declaration is shorter;
   its defining instances are defined just as overloaded functions would,
   since each of them needs to test a specific value of |testMode|.
*/


// Add to the empty mode the test commands that require that mode.
template<>
void addTestCommands<emptymode::EmptymodeTag>
  (commands::CommandMode& mode, emptymode::EmptymodeTag)
{
  if (testMode == EmptyMode)
    mode.add("test",test_f);

  mode.add("extract-graph",extract_graph_f);
  mode.add("extract-cells",extract_cells_f);
}


// Add to the main mode the test commands that require that mode.
template<>
void addTestCommands<mainmode::MainmodeTag>
  (commands::CommandMode& mode, mainmode::MainmodeTag)
{
  if (testMode == MainMode)
    mode.add("test",test_f);

  // add additional commands here :

  mode.add("cmatrix",cmatrix_f);
  mode.add("coroots_rootbasis",coroots_rootbasis_f);
  mode.add("poscoroots_rootbasis",poscoroots_rootbasis_f);
  mode.add("posroots_rootbasis",posroots_rootbasis_f);
  mode.add("roots_rootbasis",roots_rootbasis_f);
  mode.add("rootdatum",rootdatum_f);
  mode.add("dualkgb",dual_kgb_f<false>);
}


// Add to the real mode the test commands that require that mode.
template<>
void addTestCommands<realmode::RealmodeTag>
  (commands::CommandMode& mode, realmode::RealmodeTag)
{
  if (testMode == RealMode)
    mode.add("test",test_f);

  mode.add("components",components_f);
  mode.add("corder",corder_f);
  mode.add("kgb",kgb_f<false>);
  mode.add("kgborder",kgborder_f);
  mode.add("block",block_f<false>);
  mode.add("dualblock",dual_block_f<false>);
  mode.add("blockd",blockd_f);
  mode.add("blockorder",blockorder_f);
  mode.add("blocku",blocku_f);
  mode.add("smallkgb",kgb_f<true>);
  mode.add("smallblock",block_f<true>);
  mode.add("smalldualkgb",dual_kgb_f<true>);
  mode.add("smalldualblock",dual_block_f<true>);
  mode.add("blockstabilizer",blockstabilizer_f);
  mode.add("dualmap",dual_map_f);
  mode.add("blockwrite",blockwrite_f);
  mode.add("klbasis",klbasis_f);
  mode.add("kllist",kllist_f);
  mode.add("primkl",primkl_f);
  mode.add("klwrite",klwrite_f);
  mode.add("wcells",wcells_f);
  mode.add("wgraph",wgraph_f);

}

// Add to the help mode the test commands that require that mode.
template<> void addTestHelp<emptymode::EmptymodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      emptymode::EmptymodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == EmptyMode) {
    mode.add("test",nohelp_h);
    insertTag(t,"test",test_tag);
  }

  mode.add("extract-graph",nohelp_h);
  insertTag(t,"extract-graph",extract_graph_tag);
  mode.add("extract-cells",nohelp_h);
  insertTag(t,"extract-cells",extract_cells_tag);

  // add additional help commands here:

  // add additional command tags here:

}


// Add to the main mode the help commands for test commands with that mode
template<> void addTestHelp<mainmode::MainmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      mainmode::MainmodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == MainMode) {
    mode.add("test",nohelp_h);
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:

  // add additional command tags here:

  mode.add("cmatrix",cmatrix_h);
  mode.add("coroots_rootbasis",nohelp_h);
  mode.add("gradings",nohelp_h);
  mode.add("poscoroots_rootbasis",nohelp_h);
  mode.add("posroots_rootbasis",nohelp_h);
  mode.add("roots_rootbasis",nohelp_h);
  mode.add("rootdatum",nohelp_h);

  // add additional command tags here :

  insertTag(t,"cmatrix",cmatrix_tag);
  insertTag(t,"coroots_rootbasis",test_tag);
  insertTag(t,"gradings",test_tag);
  insertTag(t,"poscoroots_rootbasis",test_tag);
  insertTag(t,"posroots_rootbasis",test_tag);
  insertTag(t,"roots_rootbasis",test_tag);
  insertTag(t,"rootdatum",test_tag);

}


// Add to the real mode the help commands for test commands with that mode
template<> void addTestHelp<realmode::RealmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      realmode::RealmodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == RealMode) {
    mode.add("test",nohelp_h);
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:

  mode.add("components",nohelp_h);
  mode.add("corder",nohelp_h);
  mode.add("kgborder",kgborder_h);
  mode.add("block",block_h);
  mode.add("blockd",blockd_h);
  mode.add("blockorder",blockorder_h);
  mode.add("blocku",blocku_h);
  mode.add("blockstabilizer",nohelp_h);
  mode.add("primkl",primkl_h);
  mode.add("involution",nohelp_h);
  mode.add("kgb",kgb_h);
  mode.add("klbasis",klbasis_h);
  mode.add("kllist",kllist_h);
  mode.add("klwrite",klwrite_h);
  mode.add("blockwrite",blockwrite_h);
  mode.add("wcells",wcells_h);
  mode.add("wgraph",wgraph_h);


  // add additional command tags here:
  insertTag(t,"block",block_tag);
  insertTag(t,"blockd",test_tag);
  insertTag(t,"blockorder",test_tag);
  insertTag(t,"blocku",blocku_tag);
  insertTag(t,"blockstabilizer",test_tag);
  insertTag(t,"components",test_tag);
  insertTag(t,"corder",test_tag);
  insertTag(t,"kgborder",test_tag);
  insertTag(t,"primkl",test_tag);
  insertTag(t,"involution",test_tag);
  insertTag(t,"kgb",kgb_tag);
  insertTag(t,"klbasis",klbasis_tag);
  insertTag(t,"kllist",kllist_tag);
  insertTag(t,"klwrite",klwrite_tag);
  insertTag(t,"blockwrite",blockwrite_tag);
  insertTag(t,"wcells",wcells_tag);
  insertTag(t,"wgraph",wgraph_tag);

}

} // namespace test

namespace {

void block_h()

{
  io::printFile(std::cerr,"block.help",io::MESSAGE_DIR);
}

void blockd_h()

{
  io::printFile(std::cerr,"blockd.help",io::MESSAGE_DIR);
}

void blockorder_h()

{
  io::printFile(std::cerr,"blockorder.help",io::MESSAGE_DIR);
}

void blocku_h()

{
  io::printFile(std::cerr,"blocku.help",io::MESSAGE_DIR);
}

void cmatrix_h()

{
  io::printFile(std::cerr,"cmatrix.help",io::MESSAGE_DIR);
}

void kgborder_h()

{
  io::printFile(std::cerr,"kgborder.help",io::MESSAGE_DIR);
}

void primkl_h()

{
  io::printFile(std::cerr,"primkl.help",io::MESSAGE_DIR);
}

void kgb_h()

{
  io::printFile(std::cerr,"kgb.help",io::MESSAGE_DIR);
}

void klbasis_h()

{
  io::printFile(std::cerr,"klbasis.help",io::MESSAGE_DIR);
}

void kllist_h()

{
  io::printFile(std::cerr,"kllist.help",io::MESSAGE_DIR);
}

void klwrite_h()

{
  io::printFile(std::cerr,"klwrite.help",io::MESSAGE_DIR);
}

void blockwrite_h()

{
  io::printFile(std::cerr,"blockwrite.help",io::MESSAGE_DIR);
}

void wcells_h()

{
  io::printFile(std::cerr,"wcells.help",io::MESSAGE_DIR);
}

void wgraph_h()

{
  io::printFile(std::cerr,"wgraph.help",io::MESSAGE_DIR);
}

}

/*****************************************************************************

        Chapter II -- Utility functions

******************************************************************************/

namespace {

const rootdata::RootDatum& currentRootDatum()

{
  return mainmode::currentComplexGroup().rootDatum();
}

}

/*****************************************************************************

        Chapter III -- Functions for the test commands

******************************************************************************/

namespace {

  // Main mode functions


// Print the Cartan matrix on stdout.
void cmatrix_f()
{
  latticetypes::LatticeMatrix q;

  rootdata::cartanMatrix(q,currentRootDatum());
  prettyprint::printMatrix(std::cout,q);

}


// Print information about the root datum (see testprint.cpp for details).
void rootdatum_f()
{
  try {
    ioutils::OutputFile file;
    testprint::print(file,currentRootDatum());
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}


// Print the roots in the simple root coordinates.
void roots_rootbasis_f()
{
  try {
    const rootdata::RootDatum& rd = currentRootDatum();
    ioutils::OutputFile file;

    for (rootdata::RootNbr i=0; i<rd.numRoots(); ++i)
      prettyprint::printInRootBasis(file,i,rd) << std::endl;
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Print the positive roots in the simple root coordinates.
void posroots_rootbasis_f()

{
  try {
    const rootdata::RootDatum& rd = currentRootDatum();
    ioutils::OutputFile file;

    prettyprint::printInRootBasis(file,rd.posRootSet(),rd);
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Print the coroots in the simple coroot coordinates.
void coroots_rootbasis_f()
{
  try {
    const rootdata::RootDatum rd (currentRootDatum(),tags::DualTag());
    ioutils::OutputFile file;

    for (rootdata::RootNbr i=0; i<rd.numRoots(); ++i)
      prettyprint::printInRootBasis(file,i,rd) << std::endl;
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Print the positive coroots in the simple coroot coordinates.
void poscoroots_rootbasis_f()
{
  try {
    const rootdata::RootDatum rd (currentRootDatum(),tags::DualTag());
    ioutils::OutputFile file;

    prettyprint::printInRootBasis(file,rd.posRootSet(),rd);
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

  // Real mode functions

/*
  Prints the (dual) component group of the current group. We print it out
  in terms of the canonical basis of T(2)^v
*/
void components_f()
{
  const realredgp::RealReductiveGroup& G = realmode::currentRealGroup();
  const latticetypes::ComponentList& c = G.dualComponentReps();

  if (c.size() > 0)
    std::cout << "component group is (Z/2)^" << c.size() << std::endl;
  else
    std::cout << "group is connected" << std::endl;

}


// Print the Hasse diagram of the ordering of Cartan classes.
void corder_f()
{ realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  try {
    G_R.fillCartan();

    std::cout << "hasse diagram of Cartan ordering:" << std::endl;
    realredgp_io::printCartanOrder(std::cout,G_R);
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }

}


// Print the Hasse diagram of the ordering of K orbits on G/B.
void kgborder_f()
{
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();
  G.fillCartan();

  std::cout << "kgbsize: " << G.kgbSize() << std::endl;
  ioutils::OutputFile file;

  kgb::KGB kgb(G);
  kgb.fillBruhat();
  kgb_io::printBruhatOrder(file,kgb.bruhatOrder());
}


/*
  Synopsis: prints out information about the stabilizer of a representation
  under the cross action
*/
void blockstabilizer_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  try {
    G_R.fillCartan();
    size_t cn;

    // get Cartan class; abort if unvalid
    interactive::getCartanClass(cn,G_R.cartanSet(),commands::currentLine());

    const complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    realform::RealForm drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(cn),tags::DualTag());

    ioutils::OutputFile file;
    realredgp_io::printBlockStabilizer(file,G_RI.realGroup(),cn,drf);
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }


}

// Print the kgb table (if |small|) only necessary part for one block
template<bool small>
void kgb_f()
{
  try {
    realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
    G_R.fillCartan();


    if (small)
    {
      const complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
      complexredgp::ComplexReductiveGroup dGC (G_C,tags::DualTag());
      const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
      const complexredgp_io::Interface& G_I = G_RI.complexInterface();

      // get dual real form
      realform::RealForm drf;

      interactive::getInteractive
	(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

      realredgp::RealReductiveGroup dGR(dGC, drf);
      dGR.fillCartan();

      bitmap::BitMap common=blocks::common_Cartans(G_R,dGR);

      std::cout << "relevant Cartan classes: ";
      basic_io::seqPrint(std::cout,common.begin(),common.end(),",","{","}\n");

      std::cout << "partial kgb size: " <<
	G_C.cartanClasses().KGB_size(G_R.realForm(),common) << std::endl;

      ioutils::OutputFile file;
      kgb::KGB kgb(G_R,common);
      kgb_io::printKGB(file,kgb);
    }
    else
    {
      std::cout << "kgbsize: " << G_R.kgbSize() << std::endl;
      ioutils::OutputFile file;
      kgb::KGB kgb(G_R);
      kgb_io::printKGB(file,kgb);
    }
  }
  catch(error::InputError e) {
    e("aborted");
  }
}


// Print a kgb table for a dual real form.
template<bool small>
void dual_kgb_f()
{
  try {

    if (small)
    {
      realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
      G_R.fillCartan(); // here we need only generate Cartans for real form
      const complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
      complexredgp::ComplexReductiveGroup dGC (G_C,tags::DualTag());
      const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
      const complexredgp_io::Interface& G_I = G_RI.complexInterface();

      // get dual real form
      realform::RealForm drf;

      interactive::getInteractive
	(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

      realredgp::RealReductiveGroup dGR(dGC, drf);
      dGR.fillCartan();

      bitmap::BitMap common=blocks::common_Cartans(dGR,G_R);

      std::cout << "relevant Cartan classes for dual group: ";
      basic_io::seqPrint(std::cout,common.begin(),common.end(),",","{","}\n");

      std::cout << "partial kgb size: " <<
	dGC.cartanClasses().KGB_size(drf,common) << std::endl;
      ioutils::OutputFile file;

      kgb::KGB kgb(dGR,common);
      kgb_io::printKGB(file,kgb);
    }
    else
    {
      complexredgp::ComplexReductiveGroup& G_C =
	mainmode::currentComplexGroup();
      G_C.fillCartan(); // must generate all Cartans: no real form chosen

      const complexredgp_io::Interface& G_I =
	mainmode::currentComplexInterface();
      const realform::RealFormList rfl = // get list of all dual real forms
	G_C.dualRealFormLabels(G_C.mostSplit(G_C.quasisplit()));

      realform::RealForm drf;

      interactive::getInteractive(drf,G_I,rfl,tags::DualTag());

      // the complex group must be in a variable: is non-const for real group
      complexredgp::ComplexReductiveGroup dG_C(G_C,tags::DualTag());
      realredgp::RealReductiveGroup dG(dG_C,drf);
      dG.fillCartan();

      std::cout << "dual kgbsize: " << dG.kgbSize() << std::endl;
      ioutils::OutputFile file;

      kgb::KGB kgb(dG);
      kgb_io::printKGB(file,kgb);
    }
  }
  catch(error::InputError e) {
    e("aborted");
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
}

// Print the current block
template<bool small>
void block_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  try {
    G_R.fillCartan();


    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    realform::RealForm drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    blocks::Block block(G_C,G_R.realForm(),drf,small);

    ioutils::OutputFile file;
    block_io::printBlock(file,block);

  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}

// Print the dual block of the current block
template<bool small>
void dual_block_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  try {
    G_R.fillCartan(); // must fill Cartans to get |G_R.mostSplit()|

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    realform::RealForm drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    // the complex group must be in a variable: it is non-const for block
    complexredgp::ComplexReductiveGroup dG_C(G_C,tags::DualTag());
    dG_C.fillCartan(drf);

    blocks::Block block(dG_C,drf,G_R.realForm(),small);

    ioutils::OutputFile file;
    block_io::printBlock(file,block);

  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}

// Print the current block with involutions in involution-reduced form
void blockd_f()
{
  using namespace block_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    RealForm drf;

    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());

    Block block(G_C,G_R.realForm(),drf);

    OutputFile file;
    printBlockD(file,block);

  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Print the unitary elements of the block.
void blocku_f()
{
  using namespace block_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    RealForm drf;

    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());

    Block block(G_C,G_R.realForm(),drf);

    OutputFile file;
    printBlockU(file,block);
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}


// Print the Hasse diagram for the Bruhat order on the current block
void blockorder_f()
{
  using namespace block_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }

  complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  Block block(G_C,G_R.realForm(),drf);

  std::cout << "block size: " << block.size() << std::endl;
  ioutils::OutputFile file;
  block.fillBruhat();
  kgb_io::printBruhatOrder(file,block.bruhatOrder());
}

// Print the correspondence of the current block with its dual block
void dual_map_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  try {
    G_R.fillCartan(); // must fill Cartans to get |G_R.mostSplit()|

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    realform::RealForm drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    // the complex group must be in a variable: it is non-const for block
    complexredgp::ComplexReductiveGroup dG_C(G_C,tags::DualTag());
    dG_C.fillCartan();

    blocks::Block block(G_C,G_R.realForm(),drf);
    blocks::Block dual_block(dG_C,drf,G_R.realForm());

    std::vector<blocks::BlockElt> v=blocks::dual_map(block,dual_block);

    std::ostringstream s("");
    basic_io::seqPrint(s,v.begin(),v.end(),", ","[","]\n");
    ioutils::OutputFile file;
    foldLine(file,s.str()," ");
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}

// Writes a binary file containing descent sets and ascent sets for block
void blockwrite_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  // reserve another BlockElt value
  const blocks::BlockElt noGoodAscent = blocks::UndefBlock-1;

  try {
    G_R.fillCartan();

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    realform::RealForm drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    std::ofstream block_out; // binary output files
    while (true)
      {
	std::string file_name= interactive::getFileName
	  ("File name for block output: ");
	if (file_name=="") break; // if no name given, don't open a file
	block_out.open(file_name.c_str(),
		       std::ios_base::out
		       | std::ios_base::trunc
		       | std::ios_base::binary);
	if (block_out.is_open()) break;
	std::cerr << "Failed to open file for writing, try again.\n";
      }

    blocks::Block block(G_C,G_R.realForm(),drf);

    unsigned char rank=block.rank(); // certainly fits in a byte

    std::cout << "Writing block data:\n";
    basic_io::put_int(block.size(),block_out);  // block size in 4 bytes
    block_out.put(rank);                        // rank in 1 byte

    { // output length data
      unsigned char max_length=block.length(block.size()-1);
      block_out.put(max_length);

      // basic_io::put_int(0,block_out); // obvious: no elements of length<0
      size_t l=0;
      for (blocks::BlockElt z=0; z<block.size(); ++z)
	while (block.length(z)>l)
	  {
	    basic_io::put_int(z,block_out); // record: z elements of length<=l
	    ++l;
	  }
      assert(l==max_length); // so max_length values are written

      // basic_io::put_int(block.size(),block_out);
      // also obvious: there are block.size() elements of length<=max_length
    }


    for (blocks::BlockElt y=0; y<block.size(); ++y)
      {
	bitset::RankFlags d;
	for (size_t s = 0; s < rank; ++s)
	  {
	    descents::DescentStatus::Value v = block.descentValue(s,y);
	    if (descents::DescentStatus::isDescent(v)) d.set(s);
	  }
	basic_io::put_int(d.to_ulong(),block_out); // write d as 32-bits value
      }

    for (blocks::BlockElt x=0; x<block.size(); ++x)
      {
#if VERBOSE
	std::cerr << x << '\r';
#endif
	for (size_t s = 0; s < rank; ++s)
	  {
	    descents::DescentStatus::Value v = block.descentValue(s,x);
            if (descents::DescentStatus::isDescent(v)
		or v==descents::DescentStatus::ImaginaryTypeII)
	      basic_io::put_int(noGoodAscent,block_out);
	    else if (v == descents::DescentStatus::RealNonparity)
	      basic_io::put_int(blocks::UndefBlock,block_out);
	    else if (v == descents::DescentStatus::ComplexAscent)
	      basic_io::put_int(block.cross(s,x),block_out);
	    else if (v == descents::DescentStatus::ImaginaryTypeI)
	      basic_io::put_int(block.cayley(s,x).first,block_out);
	    else assert(false);
	  }
      }
    std::cout<< "\nDone.\n";
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}

/* For each element $y$ in the block, outputs the list of non-zero K-L
   polynomials $P_{x,y}$.

   This is what is required to write down the K-L basis element $c_y$.
*/
void klbasis_f()
{
  using namespace basic_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace kl;
  using namespace kl_io;
  using namespace klsupport;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    RealForm drf;

    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());

    Block block(G_C,G_R.realForm(),drf);

    KLSupport kls(block);
    kls.fill();

    KLContext klc(kls);
    klc.fill();

    OutputFile file;
    file << "Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:"
	 << std::endl << std::endl;
    kl_io::printAllKL(file,klc);
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}


/*
  Synopsis: outputs the list of all distinct Kazhdan-Lusztig-Vogan polynomials
*/
void kllist_f()
{
  using namespace basic_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace kl;
  using namespace kl_io;
  using namespace klsupport;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    RealForm drf;

    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());

    Block block(G_C,G_R.realForm(),drf);

    KLSupport kls(block);
    kls.fill();

    KLContext klc(kls);
    klc.fill();

    OutputFile file;
    printKLList(file,klc);
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

/*!
  \brief Prints out the list of all K-L polynomials for primitive pairs.

  Explanation: x is primitive w.r.t. y, if any descent for y is also a
  descent for x, or a type II imaginary ascent. Ths means that none of
  the easy recursion formulas applies to P_{x,y}.
*/

void primkl_f()
{
  using namespace basic_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace kl;
  using namespace kl_io;
  using namespace klsupport;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    RealForm drf;

    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());

    Block block(G_C,G_R.realForm(),drf);

    KLSupport kls(block);
    kls.fill();

    KLContext klc(kls);
    klc.fill();

    OutputFile file;
    file << "Kazhdan-Lusztig-Vogan polynomials for primitive pairs:"
	 << std::endl << std::endl;
    kl_io::printPrimitiveKL(file,klc);
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Write the results of the KL computations to a pair of binary files
void klwrite_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  try {
    G_R.fillCartan();

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    realform::RealForm drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    std::ofstream matrix_out, coefficient_out; // binary output files
    {
      while (true)
	{
	  std::string file_name= interactive::getFileName
	    ("File name for matrix output: ");
	  if (file_name=="") break; // if no name given, don't open a file
	  matrix_out.open(file_name.c_str(),
			    std::ios_base::out
			  | std::ios_base::trunc
			  | std::ios_base::binary);
	  if (matrix_out.is_open()) break;
	  std::cerr << "Failed to open file for writing, try again.\n";
	}

      while (true)
	{
	  std::string file_name= interactive::getFileName
	    ("File name for coefficient output: ");
	  if (file_name=="") break; // if no name given, don't open a file
	  coefficient_out.open(file_name.c_str(),
			         std::ios_base::out
			       | std::ios_base::trunc
			       | std::ios_base::binary);
	  if (coefficient_out.is_open()) break;
	  std::cerr << "Failed to open file for writing, try again.\n";
	}
    }

    blocks::Block block(G_C,G_R.realForm(),drf);

    klsupport::KLSupport kls(block); kls.fill();

    kl::KLContext klc(kls); klc.fill();

    if (matrix_out.is_open())
      {
	std::vector<unsigned int> delta(klc.size());
	std::streamoff offset=0;
	std::cout << "Writing matrix rows:\n";
	for (blocks::BlockElt y=0; y<klc.size(); ++y)
	  {
#if VERBOSE
	    std::cerr << y << '\r';
#endif
	    std::streamoff new_offset=klc.writeKLRow(y,matrix_out);
	    delta[y]=static_cast<unsigned int>((new_offset-offset)/4);
	    offset=new_offset;
	  }

	if (true) {
	// now write the values allowing rapid location of the matrix rows
	for (blocks::BlockElt y=0; y<klc.size(); ++y)
	  basic_io::put_int(delta[y],matrix_out);

	// and finally sign file as being in new format by overwriting 4 bytes
	matrix_out.seekp(0,std::ios_base::beg);
	basic_io::put_int(filekl::magic_code,matrix_out);
	}
      }
    if (coefficient_out.is_open())
      {
	std::cout << "\nWriting all polynomial coefficients:\n";
	klc.writeKLStore(coefficient_out);
	std::cout << "Done.\n";
      }
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}


// Print the cells of the W-graph of the block.
void wcells_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  try {
    G_R.fillCartan();

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    realform::RealForm drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    blocks::Block block(G_C,G_R.realForm(),drf);

    klsupport::KLSupport kls(block); kls.fill();

    kl::KLContext klc(kls); klc.fill();

    wgraph::WGraph wg(klc.rank()); kl::wGraph(wg,klc);

    wgraph::DecomposedWGraph dg(wg);

    ioutils::OutputFile file; wgraph_io::printWDecomposition(file,dg);
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

/*
  Synopsis: outputs the W-graph corresponding to a block.
*/
void wgraph_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  try {
    G_R.fillCartan();

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    realform::RealForm drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    blocks::Block block(G_C,G_R.realForm(),drf);

    klsupport::KLSupport kls(block); kls.fill();

    kl::KLContext klc(kls); klc.fill();

    wgraph::WGraph wg(klc.rank()); kl::wGraph(wg,klc);

    ioutils::OutputFile file; wgraph_io::printWGraph(file,wg);

  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

void extract_graph_f()
{
  try {
    ioutils::InputFile block_file("block information");
    ioutils::InputFile matrix_file("matrix information");
    ioutils::InputFile polynomial_file("polynomial information");
    ioutils::OutputFile file;

    wgraph::WGraph wg=wgraph::wGraph(block_file,matrix_file,polynomial_file);
    wgraph_io::printWGraph(file,wg);
  }
  catch (error::InputError e) {
    e("aborted");
  }
}

void extract_cells_f()
{
  try {
    ioutils::InputFile block_file("block information");
    ioutils::InputFile matrix_file("matrix information");
    ioutils::InputFile polynomial_file("polynomial information");
    ioutils::OutputFile file;

    wgraph::WGraph wg=wgraph::wGraph(block_file,matrix_file,polynomial_file);
    wgraph::DecomposedWGraph dg(wg);
    wgraph_io::printWDecomposition(file,dg);
  }
  catch (error::InputError e) {
    e("aborted");
  }
}


/*
  Function invoked by the "test" command.
*/
void test_f()
{
  // put your code here, and define testMode at top of file appropriately

  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  try {
    G_R.fillCartan(); // must fill Cartans to get |G_R.mostSplit()|

    complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
    const realredgp_io::Interface& G_RI = realmode::currentRealInterface();
    const complexredgp_io::Interface& G_I = G_RI.complexInterface();

    // get dual real form
    realform::RealForm drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    // the complex group must be in a variable: it is non-const for block
    complexredgp::ComplexReductiveGroup dG_C(G_C,tags::DualTag());
    dG_C.fillCartan();

    blocks::Block block(G_C,G_R.realForm(),drf);
    blocks::Block dual_block(dG_C,drf,G_R.realForm());

    klsupport::KLSupport kls(block); kls.fill();
    kl::KLContext klc(kls); klc.fill();


    klsupport::KLSupport dual_kls(dual_block); dual_kls.fill();
    kl::KLContext dual_klc(dual_kls); dual_klc.fill();

    std::cout << ( kltest::dualityVerify(klc,dual_klc) ? "Succes" : "Failure" )
	      << std::endl;
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}

} // namespace

} // namespace atlas
