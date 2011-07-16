/*!
\file
\brief Class definitions and function declarations for the RootDatum class.
*/
/*
  This is subdatum.h

  Copyright (C) 2010 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef SUBDATUM_H  /* guard against multiple inclusions */
#define SUBDATUM_H

#include "rootdata.h" // derive from |RootSystem|
#include "prerootdata.h"
#include "weyl.h"
#include "tits.h"
#include "gradings_fwd.h"
#include "realredgp_fwd.h"
#include "kgb_fwd.h"

namespace atlas {

namespace subdatum {

/* The following data type is specific for representation theory. It is
   associated to a subsystem of the dual root datum (with positive roots
   contained in the positive parent coroots), and is derived from |RootSystem|
   at that dual side. It remains however attached to the parent root _datum_,
   contains and exports a reference to that rootdatum, and has a method to
   produce a |PreRootDatum| for the subsystem of the parent defined by simple
   coroots for the subsystem (not a full root datum, for efficientcy reasons).
 */

// A subsystem on the dual side of a given root datum
class SubSystem : public RootSystem // new system, subsytem of dual
{
  const RootDatum& rd; // parent root datum
  const WeylGroup sub_W; // Weyl group no reference: built by contructor
  RootNbrList pos_map; // map positive roots to root number in parent
  RootNbrList inv_map; // partial map back from all parent roots

  struct root_info
  { weyl::Generator simple;
    WeylWord to_simple;
    WeylWord reflection;

  root_info() : simple(~0), to_simple(), reflection() {}
  };
  std::vector<root_info> sub_root;

 public:
  SubSystem(const RootDatum& parent,
	    const RootNbrList& sub_sys // list of simple roots in subsys
           );

  static SubSystem integral // pseudo contructor for integral system
  (const RootDatum& parent, const RatWeight& gamma);

  SubSystem(const SubSystem& s) // copy contructor (used by pseudo contructor)
  : RootSystem(s) // copy base object
  , rd(s.rd) // share this one
  , sub_W(s.cartanMatrix()) // reconstruct (Weyl group cannot be copied)
  , pos_map(s.pos_map), inv_map(s.inv_map), sub_root(s.sub_root) // copy those
  {
    assert(false); // should never be actually called, but exist nonetheless
  }

  const RootDatum& parent_datum() const { return rd; }
  const WeylGroup& Weyl_group() const { return sub_W; }

  weyl::Twist twist(const WeightInvolution& theta,
		    WeylWord& ww) const;
  // output value |ww| gives |-^theta| as twisted involution for |sub|

  weyl::Twist parent_twist(const WeightInvolution& theta,
			   WeylWord& ww) const;
  // similar, but |ww| is (for subsystem) on parent side, and anchored at its
  // distinguished involution |delta|: on has |theta=parent(ww).delta|.

  PreRootDatum pre_root_datum() const; // viewed from parent side

  RootNbr parent_nr_simple(weyl::Generator s) const
  { return pos_map[s]; }

  RootNbr parent_nr(RootNbr alpha) const;

  weyl::Generator simple(weyl::Generator s) const
  { return sub_root[s].simple; } // parent simple root conjugated to |sub.s|

  const WeylWord& to_simple(weyl::Generator s) const
  { return sub_root[s].to_simple; } // parent conjugating word for |simple(s)|

  const WeylWord& reflection(weyl::Generator s) const
  { return sub_root[s].reflection; } // parent reflection corresponding to |s|

  Coweight sub_2rho() const { return rd.dual_twoRho(pos_map); }
  Weight parent_sub_2rho() const { return rd.twoRho(pos_map); }

  // untwisted action of |ww|, as matrix on parent side
  LatticeMatrix action_matrix(const WeylWord& ww) const;

  InvolutionData involution_data (const WeightInvolution& theta) const;

  Grading induced(Grading base_grading) const;

}; // |class SubSystem|


/* the class |SubDatum| is mathematically much richer than |SubSystem| (also,
   note the latter holds a reference to a parent root datum, so the
   terminology is somewhat misleading). While |SubSystem| may be regarded as a
   companion to |RootDatum|, the |SubSystem| class depends on gkmod stuff. Its
   unique constructor requires a root datum involution |theta| (coded by a KGB
   element for a real form) and an infinitesimal character. The latter defines
   a subsystem by integrality, while |theta| defines a subsystem twist and a
   word that expresses it.
 */
class SubDatum : public SubSystem
{
  WeylWord base_ww;  // we need this variable mostly in the constructor!
  WeightInvolution delta; // together with twist: what was missing
  TitsGroup Tg; // twisted Weyl group, plus stuff our base class knows
  WeylElt ini_tw; // records involution of initial |x| wrt |SubSystem|

  size_t rank() const; // forbid using this directly
 public:
  SubDatum(RealReductiveGroup& GR,
	   const RatWeight& gamma,
	   KGBElt x);

  const TitsGroup& Tits_group() const {return Tg; }
  TwistedInvolution init_twisted() const { return ini_tw; }
  const WeylWord& base_twisted_in_parent() const { return base_ww; }

  size_t semisimple_rank() const { return SubSystem::rank(); }

  WeightInvolution involution(TwistedInvolution tw) const;

}; // |class SubDatum|

} // |namespace subdatum|

} // |namespace atlas|

#endif
