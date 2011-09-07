#ifndef HKL_H
#define HKL_H

#include "blocks.h"
#include "polynomials.h"
#include <map>

namespace atlas {
namespace kl {
    typedef int hKLCoeff;
    typedef hKLCoeff hMuCoeff;
    typedef polynomials::Polynomial<hKLCoeff> hKLPol;
    typedef std::pair<std::vector<BlockElt>,std::vector<hMuCoeff> > hMuRow;
    typedef std::vector<hKLPol> hKLStore;
    typedef hKLStore::const_reference hKLPolRef;
}

namespace klsupport {

// support structure for computing the twisted KL polynomials
// for an hblock. Essentially just handles the (block dependent)
// primitivization of rows and elements as well as encapsulates
// access to the block functionality
class hKLSupport {
// interface
public:
	// constructors and destructors
	hKLSupport(hBlock& hb);
	
	// accessors
	// return a constant reference to the underlying hblock
	const hBlock& getBlock() const { return d_hblock; }

	// compute the cross action in the hblock
	BlockElt cross(size_t s, BlockElt j) const { return d_hblock.hcross(s,j); }

	// return the descent set for an element in an hblock
	const RankFlags& descentSet(BlockElt z) const { return d_descent[z]; }
	
	// get the rank, size, or length of an element
	size_t getRank() const { return d_hblock.hrank(); }
	size_t getSize() const { return d_hblock.hsize(); }
	size_t getLength(BlockElt x) const { return d_hblock.length(d_hblock.hfixed(x)); }

	// primitivization of rows or elements
	void primitivize(BitMap& pmap, const RankFlags& ds) const;
	BlockElt primitivize(BlockElt x, const RankFlags& ds) const;

	// returns the first element of the block that has length at least len
	BlockElt ll(size_t len) const { return d_ll[len]; }

// helpers
private: 
	// manipulators
	void fill();

// data
private:
	// non-owned reference
	hBlock& d_hblock;

	// descent sets
	std::vector<RankFlags> d_descent;
	std::vector<RankFlags> d_ascent;
	std::vector<BitMap> d_downset;
	
	// length less
	std::vector<BlockElt> d_ll;
};

}

namespace kl {
  
class hKLContext : public klsupport::hKLSupport {
// interface
public:
	// constructors and destructors
	hKLContext(hBlock& hb);

	// returns P_{x,y}
	hKLPolRef getPoly(BlockElt x, BlockElt y) const;

	// returns a vector containing the distinct polynomials
	const hKLStore& getPolyList() const { return d_store; }

	// compute the twisted KL polynomials
	void fill();

// helpers
protected:
	// computes the mu correction term for the given row
	void muCorrection(std::vector<hKLPol>& klvrow, const PrimitiveRow& prow, BlockElt y, size_t s);

	// stores the nonzero mu and mu_ elements for the given row in the mu tables
	void fillMuRow(BlockElt y);

	// finds a root for direct recursion
	size_t findRoot(BlockElt y);
	
	void makePrimitiveRow(PrimitiveRow& prow, BlockElt y) const;

	// these two functions together are responsible for writing recently computed
	// polynomials into the store
    void writeRow(const std::vector<hKLPol>& klvrow, const PrimitiveRow& prow, BlockElt y);
	KLIndex insertPoly(const hKLPol& p);

// data
protected:
	// indices to zero and one polynomials
	KLIndex d_one;
	KLIndex d_zero;

	// holds the primitive elements in each row
	std::vector<PrimitiveRow> d_prim;

	// holds the indices to polynomials in the store for 
	// the primitive elements in each row
	std::vector<KLRow> d_kl;

	// the polynomial store - just a vector of polynomials
	hKLStore d_store;

	// stores the nonzero mu and mu_ values for the polynomials
	// as they are computed. Speeds up the recursion since we dont have
	// to search through a bunch of zero mu values
	std::vector<hMuRow> d_mu;
	std::vector<hMuRow> d_mu_;

	// A silly replacement for a hash table - maps each polynomial
	// to a unique index in the store. Of course this means each polynomial
	// is stored twice, which is probably not ideal. Should ultimately
	// be replaced with something better
	std::map<hKLPol, KLIndex> d_pmap;
};
}}

#endif

