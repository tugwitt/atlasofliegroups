/*
  This is permutations_def.h. This file contains the definitions of the
  template functions declared in permutations.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  This file contains the definitions of the template functions declared in
  permutations.h
*/

#include <cassert>
#include "bitmap.h"

namespace atlas {

namespace permutations {

/* we leave |pull_back| with implicit instantiation, as the types to
   substitute for |T| (|DescentStatus|, |KGB_base::KGBfields|) are complicated
   and hard/impossible to specify in the file permutations.cpp, while the
   contents of that file is not visible when those types \emph{are} available.
*/

/* Permute the entries of |v| by "pulling back" through our permutation |pi|,
   so that |result[i]==v[pi[i]]| for all |i| afterwards.
*/
template<typename T>
std::vector<T> Permutation::pull_back(const std::vector<T>& v) const
{
  assert(v.size()==size());

  const Permutation& pi=*this;
  std::vector<T> result; result.reserve(v.size());

  for (unsigned long i=0; i<v.size(); ++i)
    result.push_back(v[pi[i]]);

  return result;
}

/* The method |permute| also could apply to all kinds of (local) types,
   so it is reasonable (and necessary) to instantiate it implicitly
*/

/*
  Applies our permutation |pi| to the vector |v|. In other words, we send each
  entry v[i] to the new position v[pi[i]]; this means that afterwards for all
  |i|: |new_v[pi[i]]==old_v[i]|, or equivalently $new_v[i]=old_v[pi^{-1}[i]]$.
  This notion of permuting an arbitrary sequence is a left action of $S_n$.

  This is pulling back (right multiplication) through the inverse permutation.

  We are able to perform this permutation essentially in-place, using an
  auxiliary bitmap. Pushing forwards through a cycle does however require 3
  assignments (a swap) per step, while backwards a single assignment would do
*/
template<typename T> void Permutation::permute(std::vector<T>& v) const
{
  assert(v.size()>=size());
  const Permutation& pi=*this;
  bitmap::BitMap seen(size()); // initialized empty

  for (size_t i = 0; i < size(); ++i)
    if (not seen.isMember(i))
    {
      seen.insert(i);
      if (pi[i]!=i)
      {
	T tmp=v[i]; // avoid accessing |v[i]| all the time
	for (size_t j=pi[i]; j!=i; j=pi[j]) // cycle from |pi[i]| to |i|
	{
	  std::swap(tmp,v[j]); // transpose |tmp] with other members in order
	  seen.insert(j);
	}
	v[i]=tmp;
      }
    }
}

} // |namespace permutations|

} // |namespace atlas|
