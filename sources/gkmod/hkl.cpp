#include <iterator>
#include "hkl.h"

namespace atlas {

/////////////////////////////////////////////////////////////////////////////////////////
// hKLSupport
/////////////////////////////////////////////////////////////////////////////////////////
namespace klsupport {

// constructor
hKLSupport::hKLSupport(hBlock& hb) : d_hblock(hb) {
	// fill the support structures
	fill();
}

// Identify primitive elements in a row.
// flags in pmap the elements that are 'primitive' with respect to the
// given descent set. An element x is primitive (or extremal) for ds 
// if all of the roots in ds are also descents for x.
void hKLSupport::primitivize(BitMap& pmap, const RankFlags& ds) const {
	// local data
	size_t rank = getRank();

	// compute extremal elements
	for (weyl::Generator s=0; s<rank; s++) {
		if (ds.test(s)) pmap &= d_downset[s];
	}
}

// Find the primitivization of an element x.
// If an element x is not primitive for a descent set ds
// then at least one element of ds is an ascent for x. We follow this 
// ascent (i.e. x = cross(s,x)) until x is primitive for ds.
BlockElt hKLSupport::primitivize(BlockElt x, const RankFlags& ds) const {
	RankFlags ads = d_ascent[x] & ds;
	while (ads.any()) {
		size_t s = ads.firstBit();
		x = cross(s,x); 
		ads = d_ascent[x] & ds;
	}
	return x;
}

// Fills the main data structures for this support class. The d_ll vector
// allows us to quickly scan all elements if a given length. The descent and
// ascent sets allow us to quickly primitivize elements or rows.
void hKLSupport::fill() {
	// fill the vector ll - ll[i] holds the first element of the hblock
	// that has length at least i. It is not necessarily the case
	// that all lengths between 0 and length(size-1) occurr
	size_t rank = getRank();
	size_t size = getSize();
	size_t maxlen = getLength(size-1);
	d_ll.resize(maxlen+2);
	d_ll[0] = 0;
	size_t currlen = 0;
	for (size_t x=1; x<size; x++) {
		size_t xlen = getLength(x);
		while (xlen > currlen) {
			d_ll[++currlen] = x;
		}
	}
	d_ll[currlen+1] = size;

	// fill the descent sets
	d_descent.resize(size);
	d_ascent.resize(size);
	d_downset.resize(rank);
	for (weyl::Generator s=0; s<rank; s++) {
		d_downset[s].set_capacity(size);
		for (BlockElt x=0; x<size; x++) {
			if (cross(s,x) < x) {
				// descents
				d_downset[s].insert(x);
				d_descent[x].set(s);
			}
			else {
				// ascents
				d_ascent[x].set(s);
			}
		}
	}
}

}
	
/////////////////////////////////////////////////////////////////////////////////////////
// hKLContext
/////////////////////////////////////////////////////////////////////////////////////////
namespace kl {

// constructor
hKLContext::hKLContext(hBlock& hb) : klsupport::hKLSupport(hb), d_store(2) {
	// initialize the store to contain the polynomials 0 and 1
	hKLPol one(1);
	hKLPol zero(0);
	d_one = 1;
	d_zero = 0;
	d_pmap[zero] = 0;
	d_pmap[one] = 1;
	d_store[d_one]=one;
	d_store[d_zero]=zero;
}

// lookup a polynomial from the store
hKLPolRef hKLContext::getPoly(BlockElt x, BlockElt y) const {
  const PrimitiveRow& prow = d_prim[y];
  const KLRow& klrow = d_kl[y];

  x=primitivize(x, descentSet(y));
  if (x>y) return d_store[d_zero];

  PrimitiveRow::const_iterator xptr = std::lower_bound(prow.begin(),prow.end(),x);
  if (xptr == prow.end() || *xptr != x) return d_store[d_zero];
  
  return d_store[klrow[xptr - prow.begin()]];
}

// computes the primitive elements in the given row
void hKLContext::makePrimitiveRow(PrimitiveRow& prow, BlockElt y) const {
	// local data
	BitMap pmap(getSize());
	pmap.fill(0,ll(getLength(y)));
	pmap.insert(y);                     

	// filter out those that are not primitive
	primitivize(pmap,descentSet(y));

	// copy from bitset b to list e
	prow.reserve(prow.size()+pmap.size());
	std::copy(pmap.begin(),pmap.end(),back_inserter(prow));
}

// find a type II root that reduces length for the given row and return it
// if no such root is found, return the rank of the group. In this case
// no recursion is necessary since there are no primitive elements for y
size_t hKLContext::findRoot(BlockElt y) {
	// local data
	size_t rank = getRank();
	size_t ylen = getLength(y);
	BlockElt ymax = ll(ylen-1);
	const RankFlags& droots = descentSet(y);

	// examine the roots
	for (int s=0; s<rank; s++) {
		if (droots.test(s)) {
			// descent root - see if its type II
			if (cross(s,y) < ymax) {
				// type II - we are done
				return s;
			}
		}
	}
	
	// no type II descent roots found -
	// easy case, nothing to do
	return rank;
}

// write the polynomials in klv to the store. This function should be independent of how 
// the polynomials are actually stored/mapped. Those implementation details should remain hidden
// inside the insertPoly function
void hKLContext::writeRow(const std::vector<hKLPol>& klvrow, const PrimitiveRow& prow, BlockElt y) {
	// local data
	size_t psize = prow.size();
	d_prim[y].resize(psize);
	d_kl[y].resize(psize);
	for (size_t i=0; i<psize; i++) {
		d_prim[y][i] = prow[i];
		d_kl[y][i] = insertPoly(klvrow[i]);
	}
}

// This function uses a map as a hash table, which is probably not the best approach. The map
// uses the polynomials as a key and the index in d_store as the value. Swapping this approach 
// out for a better one should only require changing this function
KLIndex hKLContext::insertPoly(const hKLPol& p) {
	// local data
	size_t ssize = d_store.size();
	std::map<hKLPol,KLIndex>::iterator it;
	
	// see if the polynomial is new
	it=d_pmap.find(p);
	if (it != d_pmap.end()) {
		// not new
		return it->second;
	}
	
	// otherwise the entry wasn't found - add it
	d_store.push_back(p);
	d_pmap[p] = ssize;
	return ssize;
}

// fills the mu and mu_ tables for the given row. The mu value for a 
// P_{x,y} can be nonzero only if x is primitive with respect to y
// or length(y) - length(x) = 1. The mu_ value for P_{x,y} can be nonzero
// only if length(x) is exactly one less than an element for which mu is nonzero (or x is primitive)
// Therefore, traversing the list of primitive elements is sufficient to determine
// all of the nonzero mu and mu_ values, provided we are careful.
void hKLContext::fillMuRow(BlockElt y) {
	// local data
	size_t rank = getRank();
	size_t psize = d_prim[y].size()-1;
	size_t ylen = getLength(y);
	
	// make sure mu_ entries are added to the 
	// tables only one time
	std::vector<bool> mu_bits(getSize());

	// start by looking at each primitive element for y
	for (size_t i=0; i<psize; i++) {
		BlockElt x = d_prim[y][i];
		size_t xlen = getLength(x);
		size_t d = (ylen-xlen-1)/2;

		// odd difference
		if ((ylen-xlen) % 2) {
			for (size_t s=0; s<rank; s++) {
				// first check for lower elements - those with even difference
				// and add them to the mu_ tables. These elements dont have to
				// be primitive
				BlockElt z = cross(s,x);
				if (mu_bits[z]) continue;
				size_t zlen = getLength(z);
				if (zlen == xlen-1) {
					const hKLPol& p = getPoly(z,y);
					if (p.degree() == d) {
						d_mu_[y].first.push_back(z);
						d_mu_[y].second.push_back(p[d]);
						mu_bits[z] = true;
					}
				}
			}
			
			// check the current element
			hKLPol& p = d_store[d_kl[y][i]];
			if (p.degree() == d) {
				d_mu[y].first.push_back(x);
				d_mu[y].second.push_back(p[d]);
			}
		}
		
		// even difference
		else {
			if (mu_bits[x]) continue;
			hKLPol& p = d_store[d_kl[y][i]];
			if (p.degree() == d) {
				d_mu_[y].first.push_back(x);
				d_mu_[y].second.push_back(p[d]);
				mu_bits[x] = true;
			}
		}
	}

	// add the extras whose length is only one or two less than y
	for (size_t s=0; s<rank; s++) {
		BlockElt x = cross(s,y);
		size_t xlen = getLength(x);
		
		// odd length
		if ((ylen - xlen) == 1) {
			const hKLPol &p = getPoly(x,y);
			if (p[0] != 0) {
				d_mu[y].first.push_back(x);
				d_mu[y].second.push_back(p[0]);
			}
		}
		
		// even length
		else if ((ylen - xlen) == 2) {
			if (mu_bits[x]) continue;
			const hKLPol &p = getPoly(x,y);
			if (p[0] != 0) {
				d_mu_[y].first.push_back(x);
				d_mu_[y].second.push_back(p[0]);
				mu_bits[x] = true;
			}
		}
	}
}

// completes the recursive step by computing the mu correction term for each primitive element
// in prow. The labelling of elements is according to the notes of Vogan/Lusztig. Here sws is
// the current row so that s is type II for w and len(sws) = l(w) + 2.
void hKLContext::muCorrection(std::vector<hKLPol>& klvrow, const PrimitiveRow& prow, BlockElt sws, size_t s) {
	// local data
	BlockElt w = cross(s,sws);
	size_t wlen = getLength(w);
	size_t psize = prow.size()-1;
	size_t musize = d_mu[w].first.size();
	
	// start with the usual mu terms.
	// correction terms exist whenever mu is nonzero
	for (size_t i=0; i<musize; i++) {
		BlockElt y = d_mu[w].first[i];
		hKLCoeff mu = d_mu[w].second[i];
		BlockElt sxy = cross(s,y);
		size_t ylen = getLength(y);
		size_t sxylen = getLength(sxy);
		size_t diff = wlen - ylen;

		// first mu correction term from the notes
		if (sxylen < ylen) {
			// compute multiplicative factor
			size_t deg = (diff+1)/2;
			hKLPol p(deg+1, -mu);
			p[deg] = -mu;
			for (size_t j=0; j<psize; j++) {
				// scale each primitive element
				BlockElt sx = prow[j];
				BlockElt x = cross(s,sx);
				size_t xlen = getLength(x);
				if (xlen > ylen) break;
				klvrow[j] += (p * getPoly(x,y));
			}
		}

		// first part of the second mu correction term
		else if (sxylen == ylen+1) {
			// compute multiplicative factor
			size_t deg = (diff+2)/2;
			hKLPol p(deg, -mu);
			for (size_t j=0; j<psize; j++) {
				// scale each primitive element
				BlockElt sx = prow[j];
				BlockElt x = cross(s,sx);
				size_t xlen = getLength(x);
				if (xlen > ylen) break;
				klvrow[j] += (p * getPoly(x,sxy));
			}
		}
	}

	// the mu_ terms show up only in the second part
	// of the second mu_ correction term
	musize = d_mu_[w].first.size();
	for (size_t i=0; i<musize; i++) {
		BlockElt y = d_mu_[w].first[i];
		hKLCoeff mu = d_mu_[w].second[i];
		BlockElt sxy = cross(s,y);
		size_t ylen = getLength(y);
		size_t sxylen = getLength(sxy);
		size_t diff = wlen - ylen;
		if (sxylen > ylen) continue;
		
		for (size_t j=0; j<psize; j++) {
			// scale each primitive element
			BlockElt sx = prow[j];
			BlockElt x = cross(s,sx);
			size_t xlen = getLength(x);
			if (xlen > ylen) break;
			const hKLPol& p_xy = getPoly(x,y);
			if (p_xy.isZero()) continue;

			size_t deg = (diff+2)/2;
			hKLPol p(deg, -mu);
			klvrow[j] += (p * p_xy);
		}
	}

	// the final term is complicated. Based on the way the mu tables are
	// stored, it is easiest to simply identify places where mu(y,z)mu(z,w)
	// are multiply by the corresponding polynomial factor. This causes there
	// to be more polynomial multiplications that is really necessary, so
	// there is probably room for improvement here. Perhaps the mu coefficients
	// could be stored in a sparse matrix that allowed for easy traversing in
	// two directions or something.
	size_t wmusize = d_mu[w].first.size();
	for (size_t i=0; i<wmusize; i++) {
		// start with places where mu(z,w) is nonzero
		BlockElt z = d_mu[w].first[i];
		size_t zlen = getLength(z);
		BlockElt sxz = cross(s,z);
		if (getLength(sxz) > zlen) continue;
		
		// now look for places where mu(y,z) is nonzero
		hKLCoeff muzw = d_mu[w].second[i];
		size_t zmusize = d_mu[z].first.size();
		for (size_t j=0; j<zmusize; j++) {
			BlockElt y = d_mu[z].first[j];
			size_t ylen = getLength(y);
			BlockElt sxy = cross(s,y);
			if (getLength(sxy) > ylen) continue;
			
			hKLCoeff muyz = d_mu[z].second[j];
			hKLCoeff muprod = muyz*muzw;
			size_t diff = wlen - ylen;
			for (size_t k=0; k<psize; k++) {
				// for each such places, apply the polynomial
				// scale factor
				BlockElt sx = prow[k];
				BlockElt x = cross(s,sx);
				size_t xlen = getLength(x);
				if (xlen > ylen) break;
				const hKLPol& p_xy = getPoly(x,y);
				if (p_xy.isZero()) continue;

				size_t deg = (diff+2)/2;
				hKLPol p(deg, muprod);
				klvrow[k] += (p * p_xy);
			}
		}
	}
}

// computes the twisted KL polynomials one row at a time by induction on length. 
// Downward induction along the rows is not necessary. For each element y, we
// look for a type II root that reduces length. If such a root is found, we then compute
// the primitive elements for y. For each such primitive element, we apply the easy recursive
// formulas first and then apply the mu corrections. Once this is complete, we store the 
// polynomials and fill the mu tables.
//
// If no type II reducing root is found, then there are no primitive elements for y
// so there is remaining work to do.
void hKLContext::fill() {
	// local data
	size_t rank = getRank();
	size_t size = getSize();
	hKLPol q(1,1);

	// compute the twisted KL polynomials
	d_prim.resize(size);
	d_kl.resize(size);
	d_mu.resize(size);
	d_mu_.resize(size);

	// base case
	d_prim[0].push_back(0);
	d_kl[0].push_back(d_one);

	// fill the lists
	size_t maxlen = getLength(size-1);
	for (size_t i=1; i<=maxlen; i++) {
		BlockElt lmax = ll(i+1);
		for (BlockElt j=ll(i); j<lmax; j++) {	
			// find the reducing root
			size_t s = findRoot(j);

			// compute the primitive elements
			PrimitiveRow prow;
			makePrimitiveRow(prow,j);
			size_t psize = prow.size();
			std::vector<hKLPol> klvrow(psize);
			if (s != rank) {
				// apply the recursive formulas
				BlockElt sxj = cross(s,j);
				for (size_t k=0; k<psize-1; k++) {
					BlockElt z = prow[k];
					BlockElt sxz = cross(s,z);
					hKLPolRef Pszsj = getPoly(sxz,sxj);
					hKLPolRef Pzsj = getPoly(z,sxj);
					if (getLength(z) - getLength(sxz) == 1) {
						// type I: (q+1) P_{s.z, s.j} + (q^2 - q) P_{z,s.j}
						klvrow[k] = q*Pszsj + Pszsj + q*q*Pzsj - q*Pzsj;
					}
					else {
						// type II: P_{s.z, s.j} + q^2 P_{z,s.j}
						klvrow[k] = Pszsj + q*q*Pzsj;
					}
				}
			}
			// last K-L polynomial is 1
			klvrow.back() = d_store[d_one];

			if (s != rank) {
				// do mu-correction
				muCorrection(klvrow,prow,j,s);
			}

			// commit the results
			writeRow(klvrow,prow,j);

			fillMuRow(j);
		}
	}
}

}}

