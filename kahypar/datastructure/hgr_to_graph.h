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

using NodeID = HypernodeID;
#define INVALID_NODE std::numeric_limits<NodeID>::max()
using EdgeID = HyperedgeID;
using EdgeWeight = long double;
using ClusterID = PartitionID;
    
struct Edge {
    NodeID targetNode;
    EdgeWeight weight;
};

struct IncidentClusterWeight {
    ClusterID clusterID;
    EdgeWeight weight;
    
    IncidentClusterWeight(ClusterID clusterID, EdgeWeight weight) 
                         : clusterID(clusterID), weight(weight) { }
};
    
using NodeIterator = std::vector<NodeID>::const_iterator;
using EdgeIterator = std::vector<Edge>::const_iterator;
using IncidentClusterWeightIterator = std::vector<IncidentClusterWeight>::const_iterator;

    
class HgrToGraph {
    
public:
    HgrToGraph(Hypergraph& hypergraph) : _N(hypergraph.initialNumNodes()+hypergraph.initialNumEdges()),
                                         _adj_array(_N+1), _nodes(_N), _edges(), 
                                         _cluster_id(_N), _incidentClusterWeight(_N,IncidentClusterWeight(0,0.0L)),
                                         _isInNeighborVector(_N) {
        std::iota(_nodes.begin(),_nodes.end(),0);
        std::iota(_cluster_id.begin(),_cluster_id.end(),0);
        construct(hypergraph);
    }
    
    HgrToGraph(std::vector<NodeID>& adj_array, std::vector<Edge>& edges) 
               : _N(adj_array.size()-1),_adj_array(adj_array), _nodes(_N), _edges(edges), 
                 _cluster_id(_N), _incidentClusterWeight(_N,IncidentClusterWeight(0,0.0L)),
                 _isInNeighborVector(_N) {
        std::iota(_nodes.begin(),_nodes.end(),0);
        std::iota(_cluster_id.begin(),_cluster_id.end(),0);
    }
    
    std::pair<NodeIterator,NodeIterator> nodes() const {
        return std::make_pair(_nodes.begin(),_nodes.end());
    }
    
    std::pair<EdgeIterator,EdgeIterator> adjacentNodes(NodeID id) const {
        ASSERT(id < numNodes(), "NodeID " << id << " doesn't exist!");
        return std::make_pair(_edges.begin()+_adj_array[id],_edges.begin()+_adj_array[id+1]);
    }
    
    //Returns an iterator to _incidentClusterWeight, which store the incident ClusterID of NodeID node and the sum
    //of all incident nodes of NodeID node within a cluster.
    std::pair<IncidentClusterWeightIterator,IncidentClusterWeightIterator> incidentClusterWeight(NodeID node) {
        
        _isInNeighborVector.clear();
        size_t idx = 0;
        
        _incidentClusterWeight[idx] = IncidentClusterWeight(clusterID(node),0.0L);
        _isInNeighborVector[clusterID(node)] = idx++;
        for(Edge e : adjacentNodes(node)) {
            NodeID id = e.targetNode;
            EdgeWeight w = e.weight;
            ClusterID c_id = clusterID(id);
            if(_isInNeighborVector.contains(c_id)) {
                size_t i = _isInNeighborVector[c_id];
                _incidentClusterWeight[i].weight += w;
            }
            else {
                _incidentClusterWeight[idx] = IncidentClusterWeight(c_id,w);
                _isInNeighborVector[c_id] = idx++;
            }
        }
        
        return std::make_pair(_incidentClusterWeight.begin(),_incidentClusterWeight.begin()+idx);
    }
    
     //Returns an iterator to _incidentClusterWeight, which store the incident ClusterID of Cluster cid and the sum
    //of all incident edge weights of ClusterID cid.
    std::pair<IncidentClusterWeightIterator,IncidentClusterWeightIterator> incidentNeighborClusters(ClusterID cid) {
        
        _isInNeighborVector.clear();
        size_t idx = 0;
        
        for(NodeID node : nodes()) {
            if(clusterID(node) == cid) {
                for(Edge e : adjacentNodes(node)) {
                    NodeID id = e.targetNode;
                    EdgeWeight w = e.weight;
                    ClusterID c_id = clusterID(id);
                    if(_isInNeighborVector.contains(c_id)) {
                        size_t i = _isInNeighborVector[c_id];
                        _incidentClusterWeight[i].weight += w;
                    }
                    else {
                        _incidentClusterWeight[idx] = IncidentClusterWeight(c_id,w);
                        _isInNeighborVector[c_id] = idx++;
                    }
                }
            }
        }
        
        for(auto cur = _incidentClusterWeight.begin(); cur < _incidentClusterWeight.begin()+idx; ++cur) {
                if(cur->clusterID == cid) {
                    cur->weight /= 2.0L; break;
                }
        }
        
        return std::make_pair(_incidentClusterWeight.begin(),_incidentClusterWeight.begin()+idx);
    }
    
    std::pair<HgrToGraph,std::vector<NodeID>> contractCluster() {
            std::vector<NodeID> cluster2Node(numNodes(),INVALID_NODE);
            std::vector<NodeID> node2contractedNode(numNodes(),INVALID_NODE);
            ClusterID new_cid = 0;
            for(NodeID node : nodes()) {
                ClusterID cid = clusterID(node);
                if(cluster2Node[cid] == INVALID_NODE) {
                    cluster2Node[cid] = new_cid++;
                }
                node2contractedNode[node] = cluster2Node[cid];
                setClusterID(node,node2contractedNode[node]);
            }
            
            std::vector<NodeID> new_adj_array(new_cid+1,0);
            std::vector<Edge> new_edges;
            for(ClusterID cid = 0; cid < new_cid; ++cid) {
                new_adj_array[cid] = new_edges.size();
                for(auto incidentClusterWeight : incidentNeighborClusters(cid)) {
                    Edge e;
                    e.targetNode = static_cast<NodeID>(incidentClusterWeight.clusterID);
                    e.weight = incidentClusterWeight.weight;
                    new_edges.push_back(e);
                }
            }
            new_adj_array[new_cid] = new_edges.size();
            
            HgrToGraph graph(new_adj_array,new_edges);
            return std::make_pair(std::move(graph),node2contractedNode);
    }
    
    size_t numNodes() const {
        return static_cast<size_t>(_N);
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
        
        std::cout << "Number Nodes: " << numNodes() << std::endl;
        std::cout << "Number Edges: " << numEdges() << std::endl;
        
        for(NodeID n : nodes()) {
            std::cout << "Node ID: " << n << ", Adj. List: ";
            for(Edge e : adjacentNodes(n)) {
                std::cout << "(" << e.targetNode << ",w=" << e.weight<<") ";
            }
            std::cout << "\n";
        }
        
    }
    
    void construct(Hypergraph& _hg) {
        NodeID sum_edges = 0;
        
        //Construct adj. array for all hypernode. Amount of edges is equal to the degree of the corresponding hypernode.
        for(HypernodeID hn : _hg.nodes()) {
            _adj_array[hn] = sum_edges;
            sum_edges += _hg.nodeDegree(hn);
        }
        
        //Construct adj. array for all hyperedges. Amount of edges is equal to the size of the corresponding hyperedge.
        for(HyperedgeID he : _hg.edges()) {
            _adj_array[mapHyperedgeToGraphNode(_hg,he)] = sum_edges;
            sum_edges += _hg.edgeSize(he);
        }
        
        _adj_array[_N] = sum_edges;
        _edges.resize(sum_edges);
        
        for(HypernodeID hn : _hg.nodes()) {
            size_t pos = 0;
            for(HyperedgeID he : _hg.incidentEdges(hn)) {
                Edge e;
                e.targetNode = mapHyperedgeToGraphNode(_hg,he);
                e.weight = static_cast<EdgeWeight>(_hg.edgeWeight(he))/static_cast<EdgeWeight>(_hg.edgeSize(he));
                _edges[_adj_array[hn] + pos++] = e;
            }
        }
        
        for(HyperedgeID he : _hg.edges()) {
           size_t pos = 0;
           for(HypernodeID hn : _hg.pins(he)) {
                Edge e;
                e.targetNode = hn;
                e.weight = static_cast<EdgeWeight>(_hg.edgeWeight(he))/static_cast<EdgeWeight>(_hg.edgeSize(he));
                _edges[_adj_array[mapHyperedgeToGraphNode(_hg,he)]+pos++] = e;
           }
        }
        
        //printGraph();
        
    }
    
    inline NodeID mapHyperedgeToGraphNode(const Hypergraph& _hg, const NodeID he) const {
        return _hg.initialNumNodes() + he;
    }
    
    NodeID _N;
    std::vector<NodeID> _adj_array;
    std::vector<NodeID> _nodes;
    std::vector<Edge> _edges;
    std::vector<ClusterID> _cluster_id;
    std::vector<IncidentClusterWeight> _incidentClusterWeight;
    
    SparseMap<ClusterID,size_t> _isInNeighborVector;
    
};
        
        
}  // namespace ds
}  // namespace kahypar

