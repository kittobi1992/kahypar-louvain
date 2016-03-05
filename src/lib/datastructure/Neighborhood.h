/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_LIB_DATASTRUCTURE_NEIGHBORHOOD_H_
#define SRC_LIB_DATASTRUCTURE_NEIGHBORHOOD_H_

#include <vector>
#include <algorithm>

#include "lib/definitions.h"
#include "lib/datastructure/SparseSet.h"

using defs::Hypergraph;
using defs::HypernodeID;

namespace datastructure {

class Neighborhood {
 public:
  explicit Neighborhood(const Hypergraph& hypergraph) :
      _neighbors(hypergraph.initialNumNodes()){
    SparseSet<HypernodeID> neighbors(hypergraph.initialNumNodes());

    for (const auto hn : hypergraph.nodes()) {
      for (const auto he : hypergraph.incidentEdges(hn)) {
        for (const auto pin : hypergraph.pins(he)) {
          neighbors.add(pin);
        }
      }

      _neighbors[hn].swap(neighbors.elements());
      std::sort(_neighbors[hn].begin(), _neighbors[hn].end());

      neighbors.clear();
    }
  }

  const std::vector<HypernodeID> & of(const HypernodeID hn) const {
    return _neighbors[hn];
  }

  HypernodeID sizeOfIntersection(const HypernodeID u, const HypernodeID v) const {
    std::vector<HypernodeID> intersection;
    std::set_intersection(_neighbors[u].begin(),
                   _neighbors[u].end(),
                   _neighbors[v].begin(),
                   _neighbors[v].end(), std::back_inserter(intersection));
    return intersection.size();
  }

  HypernodeID sizeOfUnion(const HypernodeID u, const HypernodeID v) const {
    std::vector<HypernodeID> set_union;
    std::set_union(_neighbors[u].begin(),
                   _neighbors[u].end(),
                   _neighbors[v].begin(),
                   _neighbors[v].end(), std::back_inserter(set_union));
    return set_union.size();
  }

  double jaccardIndex(const HypernodeID u, const HypernodeID v) const {
    return static_cast<double>(sizeOfIntersection(u,v)) / sizeOfUnion(u,v);
  }

  void mergeNeighbors(const HypernodeID u, const HypernodeID v) {
    std::vector<HypernodeID> set_union;
    std::set_union(_neighbors[u].begin(),
                   _neighbors[u].end(),
                   _neighbors[v].begin(),
                   _neighbors[v].end(), std::back_inserter(set_union));
    _neighbors[u].swap(set_union);
  }

  std::vector<std::vector<HypernodeID>> _neighbors;

};

}  // namespace datastructure
#endif  // SRC_LIB_DATASTRUCTURE_NEIGHBORHOOD_H_
