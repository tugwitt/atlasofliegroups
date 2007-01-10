/*!
\file
\brief
Class definitions and function declarations for the class KLContext.
*/
/*
  This is kl.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef KL_H  /* guard against multiple inclusions */
#define KL_H

#include <limits>
#include <pthread.h>
#include <iostream>

#include "kl_fwd.h"

#include "blocks_fwd.h"

#include "bitset.h"
#include "klsupport.h"

#include "hashtable-stats.h-thread"
#include "polynomials.h"
#include "prettyprint.h"
#include "wgraph.h"

namespace atlas {

  /* Chapter 0
     Preliminary definitions, needed before KLContext can be defined
  */

// this one is really only to facilitate tracing, using "hashtable-loud.h"
std::ostream& operator<<  (std::ostream& out, const kl::KLPol& p);

namespace kl {

  /*!
\brief Wrapper for passing member function fillKLRow to pthread_create.
  */
  extern "C" void *ThreadStartup(void *);

  extern size_t NThreads;

// storage for KL polynomials

/* The idea is to store only sequences of coefficients in a huge vector
   |pool|, and to provide minimal overhead (the |index| vector) to be able to
   access individual polynomials. Polynomials are extracted in the form of a
   |polynomials::PolRef<KLCoeff>| object, which class replaces the type |const
   KLPol&| that is impossible to produce without having an actual |KLPol|
   value in memory. This explains our typedef of |const_reference|. Insertion
   is done by the method |push_back|, and extraction by |operator[]|,
   mimicking the behaviour of |std::vector<KLPol>| except for the return type
   from |operator[]|. Other methods are |size|, needed by the hash table code
   to know the sequence number to assign to a new polynomial sent to the
   storage, and |swap| for the |swap| method of hash tables. The methods
   |mem_size| and |mem_capacity| give the numbers of bytes used respectively
   reserved in storage, for statistical convenience.

*/

class KLPool
  {
    // a packed structure of size 1
    struct packed_byte
    {
      unsigned char b;

      unsigned int degree() const    { return b&0x1F; } // lower 5 bits
      unsigned int valuation() const { return b>>5; }   // upper 3 bits

      packed_byte(unsigned int d, unsigned int v) : b((d&0x1F) | (v<<5)) {}
      packed_byte() : b(0) {} // to initialise array members
    };

    // some parameters for storage layout
    static const int int_bits  =std::numeric_limits<unsigned int>::digits;
    static const size_t low_mask; // bit mask for lower order unsigned int

    static const unsigned int deg_limit=32; // (hard) degree limit
    static const unsigned int val_limit=8;  // (soft) valuation limit
    static const int group_bits=4; // log_2 of nr of indices packed in a group
    static const int group_size= 1<<group_bits; // nr of indices in a group
    static const unsigned int group_mask=group_size-1; // matching bit mask

/* The following kludge is necessary to repair the broken operation of bit
   shift operations with shift amount equal to the integer's width. Apparently
   only the final bits of the shift amount are used, so that x>>int_bits gives
   back x if it has exactly int_bits bits rather than setting the result to 0.
   With these definitions the right thing is done one 64-bit machines, while
   on 32-bit machines, although the high order byte we want to store is not
   present, at least it is set to 0 so that it does not perturb the (low
   order) value either, and the code can be tested as long as we do not run
   out of memory for coefficient storage (on a 32-bit machine, this will
   certainly happen before we have accumulated 2^32 coefficents).

   Unfortunately g++ warns about "shift amount >= word size" when compiling
   for 32 bits machines, even though the instructions doing that are never
   reached there (and hopefully optimised away). We have some advice for g++
   implementers here: if the compiler detects this too large shift amount,
   necessarily constant, why still produce the bitwise shift instruction,
   whose effect (returning x for x>>32 or x<<32) is very unlikely to be what
   the programmer intended, instead of generating a constant 0 which is much
   more logical. Yes, this would means that the effect of shifting over the
   same large amount would depend on whether that amount is a compile-time
   constant (and triggers a warning) or not, but the fact that the behaviour
   is officially undefined anyway gives implementers the licence to do this.

   On a 32-bit machine, set_high_order will only be called with x==0, so one
   could just write size_t(x)<<int_bits in that function to avoid the test,
   after verifying that 0<<32 == 0
*/

    static inline unsigned int high_order_int(size_t x)
      { return sizeof(x)>sizeof(unsigned int) ? x>>int_bits : 0; }
    static inline size_t set_high_order(unsigned int x)
      { return sizeof(size_t)>sizeof(unsigned int) ? size_t(x)<<int_bits : 0; }


/* The following structure collects information about a group of |group_size|
   consecutive polynomials. The address of the very first coefficient is
   recorded in 32+5=37 bits, and for all but the last polynomial the degree
   |deg<32| and a valuation |val<8| are stored in a |packed_byte|, which
   together detemine the number |1+deg-val| of stored coefficients. The
   valuation of the last polynomial of the group will be stored in the
   |pool_index_high| field of the next |IndexType| structure, and its number
   of coefficients is implicitly determined by the number of coefficients
   remaining between the first coefficient of the current group and that of
   the next one.
*/

    struct IndexType
    {
      // data members
      unsigned int pool_index_low; // lower order 32 bits of index into pool
      packed_byte pool_index_high; // high order 5 bits of above, +3 bits val
      packed_byte deg_val[group_size-1]; // degree/valuation of polynomials

      IndexType(size_t i, unsigned int v)
	: pool_index_low(i&low_mask)
        , pool_index_high(high_order_int(i),v) {}
    };

/* In a departure from the non-threaded model, we shall not directly store
   vectors of coefficients and |IndexType| blocks, but rather vectors of such
   vectors. The point is that we wish to avoid automatic reallocation at any
   time, so instead we allocate in fixed-capacity blocks, and access them via
   a fixed-capacity vector of such blocks; the only allocation done is by
   calling the reserve method for block before the first use is made of it. As
   a consequence a maximal amout of available storage is fixed a priori, but
   we shall choose a generous (but hard) limit.
*/

    std::vector<std::vector<KLCoeff> > pool;
    std::vector<std::vector<IndexType> > index;
    unsigned int pool_block_bits;  // number of bits for index within block
    unsigned int index_block_bits; // number of bits for index within block
    size_t pool_block_mask;  // mask with final pool_block_bits bits set
    size_t index_block_mask; // mask with final index_block_bits bits set

    unsigned int last_group_size; // nr of bytes of last index struct in use

    size_t savings; // gather statistics about savings by using valuations

  public:
    typedef polynomials::PolRef<KLCoeff> const_reference; // used by HashTable

    // constructor and destructor
    KLPool(size_t nr_pool_blocks, unsigned int p_block_bits,
	   size_t nr_index_blocks, unsigned int i_block_bits);
    ~KLPool(); // deletes pointers in pool and index, may print statistics

    // accessors
    const_reference operator[] (KLIndex i) const; // select polynomial by nr

    size_t size() const       // number of entries
      { return
	  ((index.size()-1<<index_block_bits) // for full blocks
	   +index.back().size()-1   // for full index groups in last block
	   <<group_bits             // multiply both these by group size
          )
	  + last_group_size; // add part of last index group used
      }

    // accessors for statistics
    size_t pool_size() const;     // number of coefficient bytes;
    size_t mem_size() const;      // net memory footprint
    size_t mem_capacity() const;  // gross memory footprint

    // manipulators
    void push_back(const KLPol&);

    void swap(KLPool& other)
      {
	pool.swap(other.pool);
	index.swap(other.index);
	std::swap(pool_block_bits,other.pool_block_bits);
	std::swap(pool_block_mask,other.pool_block_mask);
	std::swap(index_block_bits,other.index_block_bits);
	std::swap(index_block_mask,other.index_block_mask);
	std::swap(last_group_size,other.last_group_size);
	std::swap(savings, other.savings);
      }
  }; // class KLPool



typedef KLPool KLStore;

typedef KLStore::const_reference KLPolRef;

typedef KLIndex KLPtr;

typedef std::vector<KLPtr> KLRow;



} // namespace kl

/******** function declarations *********************************************/

namespace kl {

  void wGraph(wgraph::WGraph&, const KLContext&);

}

/******** type definitions **************************************************/

/* Namely: the definition of KLContext itself */


namespace kl {

 using blocks::BlockElt;

  /*!
\brief Calculates and stores the Kazhdan-Lusztig polynomials for a
block of representations of G.
  */
class KLContext {

 protected:  // permit access of our Helper class to the data members

  /*!
\brief Records whether the KL polynomials for the block have all been computed.
  */
  enum State { KLFilled, NumStates };

  /*!
\brief Bit 0 flags whether the KL polynomials have
all been computed.
  */
  bitset::BitSet<NumStates> d_state;

  /*!
\brief Pointer to the KLSupport class for this block.
  */
  klsupport::KLSupport* d_support;   // non-owned pointer

  /*!
\brief Entry d_prim[y] is a list of the elements x_i that are primitive
with respect to y and have P_{y,x_i} not zero.
  */
  std::vector<klsupport::PrimitiveRow> d_prim;

  /*!
\brief Entry d_kl[y] is a list of pointers to the polynomials
P_{y,x_i}, numbered as in the list d_prim[y].
  */
  std::vector<KLRow> d_kl;           // list of polynomial pointers

  /*!
\brief Entry d_mu[y] is a list of MuData, which are pairs (x, top
degree coefficient of P_{y,x}).
  */
  std::vector<MuRow> d_mu;           // list of mu-coefficients

  /*!
\brief Set of KL polynomials.
  */
  KLStore d_store;           // the distinct actual polynomials
  /*!
\brief Pointer to the polynomial 0.
  */
  KLPtr d_zero;
  /*!
\brief Pointer to the polynomial 1.
  */
  KLPtr d_one;

public:

// constructors and destructors
  KLContext(klsupport::KLSupport&); // initial base object with dummy pool

  // there is no point in making the destructor virtual
  ~KLContext() {}

// copy, assignment and swap
  KLContext(const KLContext&);  // ordinary copy constructor, not really used
  KLContext(const KLContext&    // augmented copy constructor, called by Helper
	    ,size_t nr_pool_blocks, size_t p_block_bits
	    ,size_t nr_index_blocks, size_t i_block_bits); // specify pool size

  KLContext& operator= (const KLContext&);

  void swap(KLContext&);

// accessors
  const blocks::Block& block() const {
    return d_support->block();
  }

  // the following two were moved here from the Helper class
  void makeExtremalRow(klsupport::PrimitiveRow& e, BlockElt y) const;

  void makePrimitiveRow(klsupport::PrimitiveRow& e, BlockElt y) const;

  /*!
\brief List of the elements x_i that are primitive with respect to y and have
 P_{y,x_i} NOT ZERO. This method is somewhat of a misnomer
  */
  const klsupport::PrimitiveRow& primitiveRow(BlockElt y) const {
    return d_prim[y];
  }

  const bitset::RankFlags& descentSet(BlockElt y) const {
    return d_support->descentSet(y);
  }

  bool isZero(const KLPtr p) const {
    return p == d_zero;
  }

  /*!
\brief The Kazhdan-Lusztig-Vogan polynomial P_{x,y}

Note that it is returned by value, since we want to allow for the possibility
of compressed storage, in chich case we cannot return an uncompressed
polynomial by reference, but we can return it by value
  */
  KLPolRef klPol(BlockElt x, BlockElt y) const;

  /*!
\brief Returns the list of pointers to the non-zero KL polynomials
P_{y,x_i} (with x_i primitive with respect to y).
  */
  const KLRow& klRow(BlockElt y) const {
    return d_kl[y];
  }

  /*!
\brief Length of y as a block element.
  */
  size_t length(BlockElt y) const {
    return d_support->length(y);
  }

  MuCoeff mu(BlockElt x, BlockElt y) const;

  /*!
\brief List of MuData, which are pairs (x, top degree coefficient of
P_{y,x}).
  */
  const MuRow& muRow(BlockElt y) const {
    return d_mu[y];
  }

  /*!
\brief Returns the set of all non-zero KL polynomials for the block.
  */
  const KLStore& polStore() const {
    return d_store;
  }

  /*!
\brief Rank of the group.
  */
  const size_t rank() const {
    return d_support->rank();
  }

  /*!
\brief Size of the block.
  */
  const size_t size() const {
    return d_kl.size();
  }

// accessors for perfoming output

  // get map of primitive elements for row y with nonzero KL polynomial
  bitmap::BitMap primMap (BlockElt y) const;

  void writeKLRow (BlockElt y, std::ostream& out) const; // write out d_kl[y]

  void writeKLStore (std::ostream& out) const;   // write out d_store

// manipulators

  // this method used to be virtual, but that seems completely silly. MvL
  void fill();


};

}

}

#endif
