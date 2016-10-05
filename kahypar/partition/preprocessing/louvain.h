/***************************************************************************
 *  Copyright (C) 2015 Tobias Heuer <tobias.heuer@gmx.net>
 *  Copyright (C) 2015-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#pragma once

#include <algorithm>
#include <limits>
#include <vector>

#include "kahypar/macros.h"
#include "kahypar/definitions.h"
#include "kahypar/datastructure/graph.h"
#include "kahypar/meta/mandatory.h"

namespace kahypar {

using ds::Graph;

template<class QualityMeasure = Mandatory>
class Louvain {
    
public:
    
    Louvain(const Hypergraph& hypergraph) : graph(hypergraph), quality(graph) { }

private:
    Graph graph;
    QualityMeasure quality;
    
};

}  // namespace kahypar

