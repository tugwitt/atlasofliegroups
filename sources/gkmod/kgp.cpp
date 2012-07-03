#include "kgp.h"

#include "constants.h"
#include "gradings.h"
#include "realredgp.h"
#include "ioutils.h"

#include <vector>
#include <bitset>
#include <algorithm>
#include <iomanip>

namespace atlas {
  namespace kgb {

    // Methods of |KGPOrbit|

std::ostream& KGPOrbit::print(std::ostream& strm) const
{
  size_t size = data.size();

  // print the orbit elements
  strm << "{";
  for (size_t i=0; i<size-1; i++)
    strm << data[i] << ",";
  strm << data[size-1] << "}";

  return strm;
}

    ///////////////////////////////////////////////////////////////////////////////////
    // KGP functions
    ///////////////////////////////////////////////////////////////////////////////////

    // Constructor
    KGP::KGP(realredgp::RealReductiveGroup& G_R, const unsigned int generators) : kgb(G_R.kgb()), kgborder(G_R.Bruhat_KGB()), kgptable(kgb.size(),(KGPElt)-1), bruhat(NULL), msize(0) {
      // local data
      size_t rank = kgb.rank();
      size_t kgbsize = kgb.size();
      std::bitset<constants::RANK_MAX> roots((int)generators);
      KGPElt unassigned = (KGPElt)-1;

      // determine the KGP orbit for each KGB orbit
      size_t count = 0;
      std::queue<KGBElt> q;
      for (KGBElt i=0; i<kgbsize; i++) {
	// see if we found a new orbit
	if (kgptable[i] == unassigned) {
	  // new orbit
	  kgptable[i] = count++;
	  q.push(i);
	}

	// for each element of the queue, explore each of its edges
	while (!q.empty()) {
	  // get the orbit
	  KGBElt kgbelt = q.front();

	  // for each simple root, determine the other elements
	  // look for other elements in the orbit
	  for (size_t j=0; j<rank; j++) {
	    if (roots[j]) {
	      // the root is in the parabolic, check the cross action
	      KGBElt ca = kgb.cross(j,kgbelt);

	      // see if we found a new element
	      if (kgptable[ca] == unassigned) {
		// new element - add it to the queue
		kgptable[ca] = kgptable[kgbelt];
		q.push(ca);
	      }

	      // if the root is noncompact, also check the Cayley transform
	      gradings::Status::Value rt = kgb.status(j,kgbelt);
	      if (rt == gradings::Status::ImaginaryNoncompact) {
		// see if the cayley transform gives a new element
		KGBElt ct = kgb.cayley(j,kgbelt);
		if (kgptable[ct] == unassigned) {
		  // new element - add it to the queue
		  kgptable[ct] = kgptable[kgbelt];
		  q.push(ct);
		}
	      }
	    }
	  }

	  // at this point we have processed the orbit
	  // so remove it from the queue
	  q.pop();
	}
      }

      // allocate enough memory to hold the orbits
      data.resize(count);

      // fill the orbits
      for (KGBElt i=0; i<kgbsize; i++) {
	data[kgptable[i]].data.push_back(i);
	if (data[kgptable[i]].size() > msize) {
	  msize = data[kgptable[i]].size();
	}
      }

      // sort by dimension - i.e. length of the open orbit
      std::sort(data.begin(), data.end());

      // the above sort undoubtedly messed up the mapping
      // so we need to fix it
      for (KGPElt i=0; i<count; i++) {
	KGPOrbit kgporbit = data[i];
	size_t osize = kgporbit.size();
	for (size_t j=0; j<osize; j++) {
	  kgptable[kgporbit.data[j]] = i;
	}
      }
    }

    // fill closure function
    void KGP::fillClosure() {
      if (bruhat != NULL) return;

      // local data
      size_t kgpsize = data.size();

      // build the Hasse diagram
      // use a bit vector to keep track of closure relations
      std::vector<set::EltList> hasse(kgpsize);
      std::vector<bool> closure(kgpsize,0);
      for (KGPElt i=0; i<kgpsize; i++) {
	// for each kgb orbit in this kgp orbit, examine closure edges
	// of degree one
	KGPOrbit kgporbit = data[i];
	size_t osize = kgporbit.size();
	size_t minelt = kgpsize;
	for (size_t j=0; j<osize; j++) {
	  // get closure edges
	  const set::EltList& clist = kgborder.hasse(kgporbit.data[j]);
	  size_t lsize = clist.size();

	  // fill the kgp closure list
	  for (size_t k=0; k<lsize; k++) {
	    KGPElt currorbit = kgptable[clist[k]];
	    closure[currorbit]=1;
	    if (currorbit < minelt) minelt = currorbit;
	  }
	}

	// at this point, closure contains a generating set of edges for
	// the closure relation. We now reduce this set to a minimal
	// generating set
	std::queue<KGPElt> q;
	for (KGPElt j=i; j-->minelt;) {
	  if (closure[j]==1) {
	    hasse[i].push_back(j);
	    q.push(j);
	    reduce(q, closure, hasse, minelt);
	  }
	}

	closure[i]=0;
	std::sort(hasse[i].begin(), hasse[i].end());
      }

      // store the hasse diagram
      bruhat = new bruhat::BruhatOrder(hasse);
    }

    // helper function - removes redundant edges from a closure relation
    void KGP::reduce(std::queue<KGPElt>& q, std::vector<bool>& closure, std::vector<set::EltList>& hasse, KGPElt minelt) {
      // while the queue is not empty, recursively remove edges
      while(!q.empty()) {
	// get the next element
	KGPElt currelt = q.front();

	// remove it from the list
	closure[currelt]=0;

	// walk the list of lower edges
	set::EltList& clist = hasse[currelt];
	size_t lsize = clist.size();
	for (size_t i=0; i<lsize; i++) {
	  if (clist[i] >= minelt) {
	    q.push(clist[i]);
	  }
	}

	// remove the element
	q.pop();
      }
    }

    // print functions
    std::ostream& KGP::print(std::ostream& strm) const {
      // local data
      size_t kgbsize = kgb.size();
      size_t kgpsize = data.size();
      size_t kgbwidth = ioutils::digits(kgbsize-1,10);
      size_t kgpwidth = ioutils::digits(kgpsize-1,10);
      size_t lwidth = ioutils::digits(kgb.length(kgbsize-1)-1,10);
      size_t cwidth = ioutils::digits(msize,10);

      // print orbits
      for (size_t i=0; i<kgpsize; i++) {
	strm << std::setw(kgpwidth) << i << ":[" << std::setw(kgbwidth) << data[i].open() << "*] ";
	strm << std::setw(lwidth) << kgb.length(data[i].open()) << " " << std::setw(cwidth) << data[i].size() << " ";
	data[i].print(strm);
	strm << std::endl;
      }

      return strm;
    }

    std::ostream& KGP::printClosure(std::ostream& strm) const {
      // local data
      size_t kgpsize = data.size();
      size_t kgpwidth = ioutils::digits(kgpsize-1,10);

      for (KGPElt i=0; i<kgpsize; i++) {
	// print the orbit
	strm << std::setw(kgpwidth) << i << ": ";

	// print the list
	const set::EltList& clist = bruhat->hasse(i);
	size_t lsize = clist.size();

	for (size_t j=0; j<lsize; j++) {
	  strm << clist[j];
	  if (j!=lsize-1) strm << ",";
	}
	strm << std::endl;
      }

      return strm;
    }

    // make a '.dot' file that can be processed by the 'dot' program
    // see www.graphviz.org for more info
    void KGP::makeDotFile(std::ostream& strm) {
      // make sure the closure order has been computed
      fillClosure();

      // local data
      size_t kgpsize = data.size();

      // write header
      strm << "digraph G {" << std::endl << "ratio=\"1.5\"" << std::endl << "size=\"7.5,10.0\"" << std::endl;

      // create vertices
      for (size_t i=0; i<kgpsize; i++) {
	// create the vertex
	strm << "v" << i << std::endl;
      }

      // add edges
      for (size_t i=0; i<kgpsize; i++) {
	const set::EltList& clist = bruhat->hasse(i);
	size_t clsize = clist.size();
	for (size_t j=0; j<clsize; j++) {
	  // add an edge in the graph
	  strm << "v" << i << " -> v" << clist[j] << "[color=gray] [arrowhead=none] [style=bold]" << std::endl;
	}
      }

      // write footer
      strm << "}" << std::endl;
    }

  }
}
