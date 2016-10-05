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

#include "gtest/gtest_prod.h"

#include "kahypar/macros.h"
#include "kahypar/definitions.h"
#include "kahypar/datastructure/sparse_map.h"

namespace kahypar {
namespace ds {

#define INVALID_NODE std::numeric_limits<NodeID>::max()
    
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

    
class Graph {
    
public:
    Graph(const Hypergraph& hypergraph) : _N(hypergraph.initialNumNodes()+hypergraph.initialNumEdges()),
                                          _adj_array(_N+1), _nodes(_N), _edges(), 
                                          _selfloop(_N,0.0L), _weightedDegree(_N,0.0L), _cluster_id(_N), 
                                          _incidentClusterWeight(_N,IncidentClusterWeight(0,0.0L)),
                                          _total_weight(0.0L), _posInIncidentClusterWeightVector(_N) {
        std::iota(_nodes.begin(),_nodes.end(),0);
        std::iota(_cluster_id.begin(),_cluster_id.end(),0);
        constructBipartiteGraph(hypergraph);
    }
    
    Graph(const std::vector<NodeID>& adj_array, const std::vector<Edge>& edges) 
                : _N(adj_array.size()-1),_adj_array(adj_array), _nodes(_N), _edges(edges), _selfloop(_N,0.0L),
                  _weightedDegree(_N,0.0L) , _cluster_id(_N), _incidentClusterWeight(_N,IncidentClusterWeight(0,0.0L)),
                  _total_weight(0.0L), _posInIncidentClusterWeightVector(_N) {
        std::iota(_nodes.begin(),_nodes.end(),0);
        std::iota(_cluster_id.begin(),_cluster_id.end(),0);
        
        for(NodeID node : nodes()) {
            for(Edge e : adjacentNodes(node)) {
                if(node == e.targetNode) {
                    _selfloop[node] = e.weight;
                }
                _total_weight += e.weight;
                _weightedDegree[node] += e.weight;
            }
        }
    }
    
    Graph(Graph&& other) : _N(std::move(other._N)),_adj_array(std::move(other._adj_array)), _nodes(std::move(other._nodes)),
                           _edges(std::move(other._edges)), _selfloop(std::move(other._selfloop)),
                           _weightedDegree(std::move(other._weightedDegree)),_cluster_id(std::move(other._cluster_id)), _incidentClusterWeight(std::move(other._incidentClusterWeight)),
                           _total_weight(std::move(other._total_weight)), 
                           _posInIncidentClusterWeightVector(std::move(other._posInIncidentClusterWeightVector)) { }
    
    Graph(const Graph& other): _N(other._N),_adj_array(other._adj_array), _nodes(other._nodes),
                               _edges(other._edges), _selfloop(other._selfloop),
                               _weightedDegree(other._weightedDegree),_cluster_id(other._cluster_id), _incidentClusterWeight(other._incidentClusterWeight),
                               _total_weight(other._total_weight), 
                               _posInIncidentClusterWeightVector(_N)  { }
                               
    Graph& operator=(const Graph&) = delete;
    Graph& operator=(Graph&&) = delete;
    
    
    
    std::pair<NodeIterator,NodeIterator> nodes() const {
        return std::make_pair(_nodes.begin(),_nodes.end());
    }
    
    std::pair<EdgeIterator,EdgeIterator> adjacentNodes(const NodeID node) const {
        ASSERT(node < numNodes(), "NodeID " << node << " doesn't exist!");
        return std::make_pair(_edges.cbegin()+_adj_array[node],_edges.cbegin()+_adj_array[node+1]);
    }
    
    size_t numNodes() const {
        return static_cast<size_t>(_N);
    }
    
    size_t numEdges() const {
        return _edges.size();
    }
    
    size_t degree(const NodeID node) const  {
        ASSERT(node < numNodes(), "NodeID " << node << " doesn't exist!");
        return static_cast<size_t>(_adj_array[node+1]-_adj_array[node]);
    }
    
    EdgeWeight weightedDegree(const NodeID node) const {
        ASSERT(node < numNodes(), "NodeID " << node << " doesn't exist!");
        return _weightedDegree[node];
    }
    
    EdgeWeight selfloopWeight(const NodeID node) const {
        ASSERT(node < numNodes(), "NodeID " << node << " doesn't exist!");
        return _selfloop[node];
    }
    
    EdgeWeight totalWeight() const {
        return _total_weight;
    }
    
    ClusterID clusterID(const NodeID node) const  {
        ASSERT(node < numNodes(), "NodeID " << node << " doesn't exist!");
        return _cluster_id[node];
    }
    
    void setClusterID(const NodeID node, const ClusterID c_id) {
        ASSERT(node < numNodes(), "NodeID " << node << " doesn't exist!");
        _cluster_id[node] = c_id;
    }
    
    /**
     * Creates an iterator to all incident Clusters of Node node. Iterator points to an 
     * IncidentClusterWeight-Struct which contains the incident Cluster ID and the sum of
     * the weight of all incident edges to that cluster.
     * 
     * @param node NodeID, which incident clusters should be evaluated
     * @return Iterator to all incident clusters of NodeID node
     */
    std::pair<IncidentClusterWeightIterator,IncidentClusterWeightIterator> incidentClusterWeightOfNode(const NodeID node);
    
    
    /**
     * Contracts the Graph based on the nodes ClusterIDs. Nodes with same ClusterID are contracted
     * in a single node in the contracted graph. Edges are inserted based on the sum of the weight
     * of all edges which connects two clusters. Also a mapping is created which maps the nodes of
     * the graph to the corresponding contracted nodes.
     * 
     * @return Pair which contains the contracted graph and a mapping from current to nodes to its
     * corresponding contrated nodes.
     */
    std::pair<Graph,std::vector<NodeID>> contractCluster();
    
protected:
    
    NodeID _N;
    std::vector<NodeID> _adj_array;
    std::vector<NodeID> _nodes;
    std::vector<Edge> _edges;
    std::vector<EdgeWeight> _selfloop;
    std::vector<EdgeWeight> _weightedDegree;
    std::vector<ClusterID> _cluster_id;
    std::vector<IncidentClusterWeight> _incidentClusterWeight;
    
    EdgeWeight _total_weight;
    
    SparseMap<ClusterID,size_t> _posInIncidentClusterWeightVector;
    
private:
    FRIEND_TEST(AGraph,DeterminesIncidentClusterWeightsOfAClusterCorrect);
    
    /*void printGraph() {
        
        std::cout << "Number Nodes: " << numNodes() << std::endl;
        std::cout << "Number Edges: " << numEdges() << std::endl;
        
        for(NodeID n : nodes()) {
            std::cout << "Node ID: " << n << ", Adj. List: ";
            for(Edge e : adjacentNodes(n)) {
                std::cout << "(" << e.targetNode << ",w=" << e.weight<<") ";
            }
            std::cout << "\n";
        }
        
    }*/
    
    /**
     * Creates an iterator to all incident Clusters of ClusterID cid. Iterator points to an 
     * IncidentClusterWeight-Struct which contains the incident ClusterID clusterID and the sum of
     * the weights of all incident edges from cid to clusterID.
     * 
     * @param cid ClusterID, which incident clusters should be evaluated
     * @return Iterator to all incident clusters of ClusterID cid
     */
    std::pair<IncidentClusterWeightIterator,IncidentClusterWeightIterator> incidentClusterWeightOfCluster(const ClusterID cid);
    
    
    void constructBipartiteGraph(const Hypergraph& hg) {
        NodeID sum_edges = 0;
        
        size_t N = static_cast<size_t>(hg.initialNumNodes());
        
        //Construct adj. array for all hypernode. Amount of edges is equal to the degree of the corresponding hypernode.
        for(HypernodeID hn : hg.nodes()) {
            _adj_array[hn] = sum_edges;
            sum_edges += hg.nodeDegree(hn);
        }
        
        //Construct adj. array for all hyperedges. Amount of edges is equal to the size of the corresponding hyperedge.
        for(HyperedgeID he : hg.edges()) {
            _adj_array[N + he] = sum_edges;
            sum_edges += hg.edgeSize(he);
        }
        
        _adj_array[_N] = sum_edges;
        _edges.resize(sum_edges);
        
        for(HypernodeID hn : hg.nodes()) {
            size_t pos = 0;
            for(HyperedgeID he : hg.incidentEdges(hn)) {
                Edge e;
                e.targetNode = N + he;
                e.weight = static_cast<EdgeWeight>(hg.edgeWeight(he))/
                           static_cast<EdgeWeight>(hg.edgeSize(he));
                _total_weight += e.weight;
                _weightedDegree[hn] += e.weight;
                _edges[_adj_array[hn] + pos++] = e;
            }
        }
        
        for(HyperedgeID he : hg.edges()) {
           size_t pos = 0;
           for(HypernodeID hn : hg.pins(he)) {
                Edge e;
                e.targetNode = hn;
                e.weight = static_cast<EdgeWeight>(hg.edgeWeight(he))/
                           static_cast<EdgeWeight>(hg.edgeSize(he));
                _total_weight += e.weight;
                _weightedDegree[N + he] += e.weight;
                _edges[_adj_array[N + he]+pos++] = e;
           }
        }
        
    }
    
};

std::pair<IncidentClusterWeightIterator,IncidentClusterWeightIterator> Graph::incidentClusterWeightOfNode(const NodeID node) {
    
    _posInIncidentClusterWeightVector.clear();
    size_t idx = 0;
    
    if(clusterID(node) != -1) {
        _incidentClusterWeight[idx] = IncidentClusterWeight(clusterID(node),0.0L);
        _posInIncidentClusterWeightVector[clusterID(node)] = idx++;
    }
    
    for(Edge e : adjacentNodes(node)) {
        NodeID id = e.targetNode;
        EdgeWeight w = e.weight;
        ClusterID c_id = clusterID(id);
        if(c_id != -1) {
            if(_posInIncidentClusterWeightVector.contains(c_id)) {
                size_t i = _posInIncidentClusterWeightVector[c_id];
                _incidentClusterWeight[i].weight += w;
            }
            else {
                _incidentClusterWeight[idx] = IncidentClusterWeight(c_id,w);
                _posInIncidentClusterWeightVector[c_id] = idx++;
            }
        }
    }
    
    return std::make_pair(_incidentClusterWeight.begin(),_incidentClusterWeight.begin()+idx);
}

std::pair<IncidentClusterWeightIterator,IncidentClusterWeightIterator> Graph::incidentClusterWeightOfCluster(const ClusterID cid) {
    
    _posInIncidentClusterWeightVector.clear();
    size_t idx = 0;
    
    for(NodeID node : nodes()) {
        if(clusterID(node) == cid) {
            for(Edge e : adjacentNodes(node)) {
                NodeID id = e.targetNode;
                EdgeWeight w = e.weight;
                ClusterID c_id = clusterID(id);
                if(_posInIncidentClusterWeightVector.contains(c_id)) {
                    size_t i = _posInIncidentClusterWeightVector[c_id];
                    _incidentClusterWeight[i].weight += w;
                }
                else {
                    _incidentClusterWeight[idx] = IncidentClusterWeight(c_id,w);
                    _posInIncidentClusterWeightVector[c_id] = idx++;
                }
            }
        }
    }
    
    return std::make_pair(_incidentClusterWeight.begin(),_incidentClusterWeight.begin()+idx);
}

std::pair<Graph,std::vector<NodeID>> Graph::contractCluster() {
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
        for(auto incidentClusterWeight : incidentClusterWeightOfCluster(cid)) {
            Edge e;
            e.targetNode = static_cast<NodeID>(incidentClusterWeight.clusterID);
            e.weight = incidentClusterWeight.weight;
            new_edges.push_back(e);
        }
    }
    new_adj_array[new_cid] = new_edges.size();
    
    Graph graph(new_adj_array,new_edges);
    return std::make_pair(std::move(graph),node2contractedNode);
}
        
}  // namespace ds
}  // namespace kahypar

