/*
  This is wgraph.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef WGRAPH_H  /* guard against multiple inclusions */
#define WGRAPH_H

#include <iostream>

#include "../Atlas.h"   // must be included before any utility headers are

#include "bitset.h"	// inlines
#include "graph.h"	// containment


namespace atlas {

/******** type declarations  (see ../Atlas.h)  ***************************/


/******** function declarations *********************************************/

namespace wgraph {

std::vector<WGraph> cells(const WGraph&);

// Functions

WGraph wGraph
  ( std::ifstream& block_file
  , std::ifstream& matrix_file
  , std::ifstream& KL_file);

}

/******** type definitions **************************************************/

namespace wgraph {

/*
  The |WGraph| class provides a shell to store a graph, edge weights
  (coefficients), and descent sets. The construction of the graph is left
  entirely to the client of this class, whence write access is given to all
  fields of thi class. It might as well have left all its members be public.
*/
class WGraph
{
  size_t d_rank;
  graph::OrientedGraph d_graph;
  std::vector<WCoeffList> d_coeff;
  std::vector<RankFlags> d_descent;

 public:

// constructors and destructors
  explicit WGraph(size_t r, size_t n);

// copy, assignment and swap: nothing needed beyond defaults

// accessors
  Partition cells(graph::OrientedGraph* p = 0) const
    { return d_graph.cells(p); }

  const WCoeffList& coeffList(graph::Vertex x) const { return d_coeff[x]; }

  const RankFlags& descent(graph::Vertex x) const { return d_descent[x]; }

  const graph::EdgeList& edgeList(graph::Vertex x) const
    { return d_graph.edgeList(x); }

  const graph::OrientedGraph& graph() const { return d_graph; }

  const size_t rank() const { return d_rank; }

  size_t size() const { return d_graph.size(); }

// manipulators
  WCoeffList& coeffList(graph::Vertex x) { return d_coeff[x]; }

  RankFlags& descent(graph::Vertex x) { return d_descent[x]; }

  graph::EdgeList& edgeList(graph::Vertex x) { return d_graph.edgeList(x); }

}; // |class WGraph|

class DecomposedWGraph
{
  typedef unsigned int cell_no;

  std::vector<WGraph> d_cell; // the strong components

  std::vector<cell_no> d_part;    // assigns strong component to each BlockElt
  std::vector< std::vector<BlockElt> > d_id; // original vertex numbers

  graph::OrientedGraph d_induced; // induced graph on cells

 public:

// constructors and destructors
  explicit DecomposedWGraph(const WGraph& wg);

// copy, assignment and swap: nothing needed beyond defaults

// accessors
  size_t rank () const { return d_cell[0].rank(); } // all ranks are equal
  size_t cellCount() const { return d_cell.size(); }
  const graph::OrientedGraph& inducedGraph() const { return d_induced; }
  const wgraph::WGraph& cell (size_t c) const { return d_cell[c]; }
  const std::vector<BlockElt>& cellMembers(size_t c) const
    { return d_id[c]; }

}; // |class DecomposedWGraph|

} // |namespace wgraph|

} // |namespace atlas|

#endif
