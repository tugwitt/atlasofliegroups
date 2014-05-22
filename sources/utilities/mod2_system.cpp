/*
  This is mod2_system.cpp: Solving linear systems over the two-element field

  Copyright (C) 2014, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "mod2_system.h"
#include <cassert>
#include <iostream>

namespace atlas {

namespace mod2_system {

std::ostream& Mod2_System::equation::print(std::ostream& f) const
{
  if (lhs.empty())
    f << '0';
  else
  {
    f << 'x' << lhs[0];
    for (unsigned long i=1; i<lhs.size(); ++i)
      f << "+x" << lhs[i];
  }
  f << " = " << (rhs ? 1 : 0) << std::endl;
  return f;
}

unsigned long Mod2_System::size() const { return pivot_index.size(); }

unsigned long Mod2_System::extend(unsigned int n)
{
  unsigned long result = size();
  pivot_index.resize(result+n,no_pivot);
  return result;
}


template <typename I> // input iterator with unsigned value type
  bool Mod2_System::add (I begin, I end, unsigned int rhs)
{
  if (not consistent())
    return false; // an already inconsistent system remains so
  bitmap::BitMap lhs(size(),begin,end); // convert iterators to bitmap

  for (unsigned int i=0; i<eqn.size(); ++i)
    if (lhs.isMember(eqn[i].lhs[0]))
      rhs += eliminate(i,lhs);

  if (lhs.empty() and rhs%2==0)
    return true; // we've ended up with a trivial equation, don't add it!

  unsigned long our_index = eqn.size();
  eqn.push_back(equation(rhs));
  equation& eq = eqn.back();
  eq.lhs.assign(lhs.begin(),lhs.end());
  if (eq.lhs.empty())
    return false; // we have just added an inconsistent equation
  pivot_index[eq.lhs[0]] = our_index;
  return true;
}

// inconsistency is recorded by |eqn.back()| having empty |lhs| (and |rhs==1|)
bool Mod2_System::consistent() const
{
  return eqn.empty() or not eqn.back().lhs.empty();
}

unsigned int
Mod2_System::eliminate (unsigned long inx, bitmap::BitMap& dst) const
{
  assert (inx<eqn.size());
  const equation& eq = eqn[inx];
  assert (dst.isMember(eq.lhs[0]));
  for (std::vector<unsigned long>::const_iterator
	 it=eq.lhs.begin(); it!=eq.lhs.end(); ++it)
    dst.flip(*it);
  return eq.rhs;
}

unsigned long Mod2_System::rank() const
{
  return consistent() ? eqn.size() : ~0ul;
}

unsigned long Mod2_System::dimension() const
{
  return consistent() ? size()-eqn.size() : ~0ul;
}

// we produce a solution by setting all free unknowns to 0
bitmap::BitMap Mod2_System::a_solution() const
{
  assert(consistent()); // no point trying for inconsistent systems
  bitmap::BitMap result (size());
  for (unsigned long inx=eqn.size(); inx-->0;) // back substitution
  {
    const equation& eq=eqn[inx];
    unsigned int val = eq.rhs ? 1 : 0;
    std::vector<unsigned long>::const_iterator it=eq.lhs.begin();
    while (++it!=eq.lhs.end()) // skip the pivot unknown itself
      if (result.isMember(*it)) // track unknowns already evaluated to 1
	++val;
    result.set_mod2(*eq.lhs.begin(),val);
  }
  return result;
} // |Mod2_System::a_solution()|

// here we do back substitution too, but only to simplify the |eqn[i]|
void Mod2_System::reduce()
{
  assert(consistent()); // no point trying for inconsistent systems
  bitmap::BitMap pivot_columns (size());
  for (unsigned long inx=eqn.size(); inx-->0;) // back substitution
  {
    equation& eq=eqn[inx];
    bitmap::BitMap row(size(),eq.lhs); // convert to bitmap for convenience

    /* eliminations in |row| will not affect any pivot column except the one
       they were intended to clear; therefore one can take a snapshot here of
       those pivot in |row| that need to be cleared */
    bitmap::BitMap to_clear = row & pivot_columns;
    pivot_columns.insert(*row.begin()); // currunt unknown is henceforth pivot
    bitmap::BitMap::iterator it = to_clear.begin();
    if (not it()) // equivalent to |to_clear.empty()|, but more efficient
      continue; // avoid doing work when nothing will change

    do
    {
      assert (pivot_index[*it]!=no_pivot and pivot_index[*it]>inx);
      eq.rhs ^= eliminate(pivot_index[*it],row);
    }
    while ((++it)());
    eq.lhs.assign(row.begin(),row.end()); // replace old by new value
  }
} // |Mod2_System::reduce|

std::vector<bitmap::BitMap> Mod2_System::solution_basis()
{
  reduce();
  std::vector<bitmap::BitMap> result;
  result.reserve(dimension());
  static const unsigned long absent = ~0ul;
  std::vector<unsigned long> result_inx(size(),absent);
  for (unsigned long i=0; i<size(); ++i)
    if (pivot_index[i]==no_pivot)
    {
      result_inx[i] = result.size();
      result.push_back(bitmap::BitMap(size()));
      result.back().insert(i);
    }
  for (unsigned long i=0; i<eqn.size(); ++i)
    for (unsigned long j=1; j<eqn[i].lhs.size(); ++j)
    {
      unsigned long k = result_inx[eqn[i].lhs[j]];
      assert(k!=absent and result[k].isMember(eqn[i].lhs[j]));
      result[k].insert(eqn[i].lhs[0]);
    }

  return result;
} // |Mod2_System::solution_basis|

// instantations
typedef std::vector<unsigned long>::iterator vec_it;
  template bool Mod2_System::add(vec_it, vec_it,unsigned int);

} // |namespace mod2_system|

} // |namespace atlas|
