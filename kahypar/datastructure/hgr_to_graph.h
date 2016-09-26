/***************************************************************************
 *  Copyright (C) 2015 Tobias Heuer <tobias.heuer@gmx.net>
 *  Copyright (C) 2015-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#pragma once

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "kahypar/macros.h"
#include "kahypar/definitions.h"

namespace kahypar {
namespace ds {
     
class HgrToGraph {

using NodeID = HypernodeID;
using ClusterID = NodeID;
using IncidentIterator = std::vector<NodeID>::const_iterator;
    
public:
    HgrToGraph(Hypergraph& hypergraph) : _N(hypergraph.initialNumNodes()), _M(hypergraph.initialNumEdges()), 
                                         _hg(hypergraph), _adj_array(_N+_M+1), _nodes(_N+_M), _edges(), _cluster_id(_N+_M) {
        std::iota(_nodes.begin(),_nodes.end(),0);
        std::iota(_cluster_id.begin(),_cluster_id.end(),0);
        construct();
    }
    
    std::pair<IncidentIterator,IncidentIterator> nodes() {
        return std::make_pair(_nodes.begin(),_nodes.end());
    }
    
    std::pair<IncidentIterator,IncidentIterator> edges(NodeID id) {
        ASSERT(id < numNodes(), "NodeID " << id << " doesn't exist!");
        return std::make_pair(_edges.begin()+_adj_array[id],_edges.begin()+_adj_array[id+1]);
    }
    
    size_t numNodes() const {
        return static_cast<size_t>(_N+_M);
    }
    
    size_t numEdges() const {
        return _edges.size();
    }
    
    size_t degree(NodeID id) {
        ASSERT(id < numNodes(), "NodeID " << id << " doesn't exist!");
        return static_cast<size_t>(_adj_array[id+1]-_adj_array[id]);
    }
    
    ClusterID clusterID(NodeID id) {
        ASSERT(id < numNodes(), "NodeID " << id << " doesn't exist!");
        return _cluster_id[id];
    }
    
    void setClusterID(NodeID id, ClusterID c_id) {
        ASSERT(id < numNodes(), "NodeID " << id << " doesn't exist!");
        _cluster_id[id] = c_id;
    }
    
private:
    
    void printGraph() {
        
        std::cout << "Number Hypernodes: " << _N << std::endl;
        std::cout << "Number Hyperedges: " << _M << std::endl;
        
        for(size_t i = 0; i <= numNodes(); ++i) {
            std::cout << _adj_array[i] << (i+1 == numNodes() ? "\n" : " ");
        }
        
        for(size_t i = 0; i < numEdges(); ++i) {
            std::cout << _edges[i] << (i+1 == numEdges() ? "\n" : " ");
        }
        
        for(NodeID n : nodes()) {
            std::cout << "Node ID: " << n << ", Adj. List: ";
            for(NodeID e : edges(n)) {
                std::cout << e << " ";
            }
            std::cout << "\n";
        }
        
    }
    
    void construct() {
        NodeID sum_edges = 0;
        
        //Construct adj. array for all hypernode. Amount of edges is equal to the degree of the corresponding hypernode.
        for(HypernodeID hn : _hg.nodes()) {
            _adj_array[hn] = sum_edges;
            sum_edges += _hg.nodeDegree(hn);
        }
        
        //Construct adj. array for all hyperedges. Amount of edges is equal to the size of the corresponding hyperedge.
        for(HyperedgeID he : _hg.edges()) {
            _adj_array[mapHyperedgeToGraphNode(he)] = sum_edges;
            sum_edges += _hg.edgeSize(he);
        }
        
        _adj_array[_N+_M] = sum_edges;
        _edges.resize(sum_edges);
        
        for(HypernodeID hn : _hg.nodes()) {
            size_t pos = 0;
            for(HyperedgeID he : _hg.incidentEdges(hn)) {
                _edges[_adj_array[hn] + pos++] = mapHyperedgeToGraphNode(he);
            }
        }
        
        for(HyperedgeID he : _hg.edges()) {
           size_t pos = 0;
           for(HypernodeID hn : _hg.pins(he)) {
                _edges[_adj_array[mapHyperedgeToGraphNode(he)]+pos++] = hn;
           }
        }
        
        //printGraph();
        
    }
    
    inline NodeID mapHyperedgeToGraphNode(const NodeID he) const {
        return _N + he;
    }
    
    NodeID _N,_M;
    Hypergraph& _hg;
    std::vector<NodeID> _adj_array;
    std::vector<NodeID> _nodes;
    std::vector<NodeID> _edges;
    std::vector<ClusterID> _cluster_id;
    
};
        
        
}  // namespace ds
}  // namespace kahypar
