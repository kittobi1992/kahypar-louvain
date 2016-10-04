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
#include "kahypar/datastructure/sparse_map.h"

namespace kahypar {
namespace ds {
     
class HgrToGraph {

using NodeID = HypernodeID;
using EdgeID = HyperedgeID;
using EdgeWeight = long double;
using ClusterID = PartitionID;
using IncidentIterator = std::vector<NodeID>::const_iterator;
using NeighborClusterWeightIterator = std::vector<std::pair<ClusterID,EdgeWeight>>::const_iterator;
    
public:
    HgrToGraph(Hypergraph& hypergraph) : _N(hypergraph.initialNumNodes()), _M(hypergraph.initialNumEdges()), 
                                         _hg(hypergraph), _adj_array(_N+_M+1), _nodes(_N+_M), _edges(), 
                                         _weight(), _cluster_id(_N+_M), _incidentClusterWeight(_N+_M,std::make_pair(-1,static_cast<EdgeWeight>(0))),
                                         _isInNeighborVector(_N+_M) {
        std::iota(_nodes.begin(),_nodes.end(),0);
        std::iota(_cluster_id.begin(),_cluster_id.end(),0);
        construct();
    }
    
    std::pair<IncidentIterator,IncidentIterator> nodes() const {
        return std::make_pair(_nodes.begin(),_nodes.end());
    }
    
    std::pair<IncidentIterator,IncidentIterator> edges(NodeID id) const {
        ASSERT(id < numNodes(), "NodeID " << id << " doesn't exist!");
        return std::make_pair(_edges.begin()+_adj_array[id],_edges.begin()+_adj_array[id+1]);
    }
    
    //Returns an iterator to _incidentClusterWeight, which store the incident ClusterID of NodeID node and the sum
    //of all incident nodes of NodeID node within a cluster.
    std::pair<NeighborClusterWeightIterator,NeighborClusterWeightIterator> incidentClusterWeight(NodeID node) {
        
        _isInNeighborVector.clear();
        size_t idx = 0;
        
        _incidentClusterWeight[idx] = std::make_pair(clusterID(node),0.0L);
        _isInNeighborVector[clusterID(node)] = idx++;
        for(EdgeID e_id = _adj_array[node]; e_id < _adj_array[node] + degree(node); ++e_id) {
            NodeID id = _edges[e_id];
            EdgeWeight w = edgeWeight(e_id);
            ClusterID c_id = clusterID(id);
            if(_isInNeighborVector.contains(c_id)) {
                size_t i = _isInNeighborVector[c_id];
                _incidentClusterWeight[i].second += w;
            }
            else {
                _incidentClusterWeight[idx] = std::make_pair(c_id,w);
                _isInNeighborVector[c_id] = idx++;
            }
        }
        
        return std::make_pair(_incidentClusterWeight.begin(),_incidentClusterWeight.begin()+idx);
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
    
    EdgeWeight edgeWeight(EdgeID edge) {
        ASSERT(id < numEdgess(), "EdgeID " << edge << " doesn't exist!");
        return _weight[edge];
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
        _weight.resize(sum_edges);
        
        for(HypernodeID hn : _hg.nodes()) {
            size_t pos = 0;
            for(HyperedgeID he : _hg.incidentEdges(hn)) {
                _edges[_adj_array[hn] + pos] = mapHyperedgeToGraphNode(he);
                _weight[_adj_array[hn] + pos++] = static_cast<EdgeWeight>(_hg.edgeWeight(he))/
                                                  static_cast<EdgeWeight>(_hg.edgeSize(he));
            }
        }
        
        for(HyperedgeID he : _hg.edges()) {
           size_t pos = 0;
           for(HypernodeID hn : _hg.pins(he)) {
                _edges[_adj_array[mapHyperedgeToGraphNode(he)]+pos] = hn;
                _weight[_adj_array[mapHyperedgeToGraphNode(he)] + pos++] = static_cast<EdgeWeight>(_hg.edgeWeight(he))/
                                                                           static_cast<EdgeWeight>(_hg.edgeSize(he));
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
    std::vector<EdgeWeight> _weight;
    std::vector<ClusterID> _cluster_id;
    std::vector<std::pair<ClusterID,EdgeWeight>> _incidentClusterWeight;
    
    SparseMap<ClusterID,size_t> _isInNeighborVector;
    
};
        
        
}  // namespace ds
}  // namespace kahypar
