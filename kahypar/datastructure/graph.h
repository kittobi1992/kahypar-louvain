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
#include "kahypar/utils/randomize.h"
#include "kahypar/partition/configuration.h"

namespace kahypar {
namespace ds {

#define INVALID_NODE std::numeric_limits<NodeID>::max()
#define INVALID_CLUSTER std::numeric_limits<ClusterID>::max()

class UnionFind {
public:
    std::vector<NodeID> parent;
        
    UnionFind(size_t n) : parent(n) {
        std::iota(parent.begin(),parent.end(),0);
    }
        
    NodeID findSet(NodeID n) { // Pfadkompression
        if (parent[n] != n) parent[n] = findSet(parent[n]);
        return parent[n];
    }
        
    void linkSets(NodeID a, NodeID b) { // Union by rank.
        parent[b] = a;
    }
        
    void unionSets(NodeID a, NodeID b) { // Diese Funktion aufrufen.
        if (findSet(a) != findSet(b)) linkSets(findSet(a), findSet(b));
    }
    
    void reset() {
        std::iota(parent.begin(),parent.end(),0);   
    }
    
};    

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
    
    Graph(const Hypergraph& hypergraph, bool use_uniform_edge_weight = false) 
                                        : _N(hypergraph.currentNumNodes()+hypergraph.currentNumEdges()),
                                          _adj_array(_N+1), _nodes(_N), _shuffleNodes(_N), _edges(), 
                                          _selfloop(_N,0.0L), _weightedDegree(_N,0.0L), _cluster_id(_N), 
                                          _incidentClusterWeight(_N,IncidentClusterWeight(0,0.0L)),
                                          _total_weight(0.0L), _posInIncidentClusterWeightVector(_N),
                                          _hypernodeMapping(hypergraph.initialNumNodes()+hypergraph.initialNumEdges(),INVALID_NODE),
                                          _unionFind(_N) {
        std::iota(_nodes.begin(),_nodes.end(),0);
        std::iota(_shuffleNodes.begin(),_shuffleNodes.end(),0);
        std::iota(_cluster_id.begin(),_cluster_id.end(),0);
        constructBipartiteGraph(hypergraph,use_uniform_edge_weight);
    }
    
    Graph(const std::vector<NodeID>& adj_array, const std::vector<Edge>& edges) 
                : _N(adj_array.size()-1),_adj_array(adj_array), _nodes(_N), _shuffleNodes(_N), _edges(edges), _selfloop(_N,0.0L),
                  _weightedDegree(_N,0.0L) , _cluster_id(_N), _incidentClusterWeight(_N,IncidentClusterWeight(0,0.0L)),
                  _total_weight(0.0L), _posInIncidentClusterWeightVector(_N), _hypernodeMapping(_N,INVALID_NODE), _unionFind(_N)  {
        std::iota(_nodes.begin(),_nodes.end(),0);
        std::iota(_shuffleNodes.begin(),_shuffleNodes.end(),0);
        std::iota(_cluster_id.begin(),_cluster_id.end(),0);
        std::iota(_hypernodeMapping.begin(),_hypernodeMapping.end(),0);
        
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
    
    Graph(const std::vector<NodeID>& adj_array, const std::vector<Edge>& edges, const std::vector<NodeID> hypernodeMapping, const std::vector<ClusterID> cluster_id) 
                : _N(adj_array.size()-1),_adj_array(adj_array), _nodes(_N), _shuffleNodes(_N), _edges(edges), _selfloop(_N,0.0L),
                  _weightedDegree(_N,0.0L) , _cluster_id(cluster_id), _incidentClusterWeight(_N,IncidentClusterWeight(0,0.0L)),
                  _total_weight(0.0L), _posInIncidentClusterWeightVector(_N), _hypernodeMapping(hypernodeMapping), _unionFind(_N)  {
        std::iota(_nodes.begin(),_nodes.end(),0);
        std::iota(_shuffleNodes.begin(),_shuffleNodes.end(),0);
        
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
                           _shuffleNodes(std::move(other._shuffleNodes)),_edges(std::move(other._edges)), _selfloop(std::move(other._selfloop)),
                           _weightedDegree(std::move(other._weightedDegree)),_cluster_id(std::move(other._cluster_id)), 
                           _incidentClusterWeight(_N,IncidentClusterWeight(0,0.0L)), _total_weight(std::move(other._total_weight)), 
                           _posInIncidentClusterWeightVector(std::move(other._posInIncidentClusterWeightVector)),
                           _hypernodeMapping(std::move(other._hypernodeMapping)), _unionFind(_N)  { }
    
    Graph(const Graph& other): _N(other._N),_adj_array(other._adj_array), _nodes(other._nodes), _shuffleNodes(other._shuffleNodes),
                               _edges(other._edges), _selfloop(other._selfloop),
                               _weightedDegree(other._weightedDegree),_cluster_id(other._cluster_id), _incidentClusterWeight(_N,IncidentClusterWeight(0,0.0L)),
                               _total_weight(other._total_weight), _posInIncidentClusterWeightVector(_N), _hypernodeMapping(other._hypernodeMapping), _unionFind(_N)   { }
                               
    Graph& operator=(Graph& other) {
        _N = other._N;
        _adj_array = other._adj_array;
        _nodes = other._nodes;
        _shuffleNodes = other._shuffleNodes;
        _edges = other._edges;
        _selfloop = other._selfloop;
        _weightedDegree = other._weightedDegree;
        _cluster_id = other._cluster_id;
        _incidentClusterWeight.assign(_N,IncidentClusterWeight(0,0.0L));
        _total_weight = other._total_weight;
        _posInIncidentClusterWeightVector = other._posInIncidentClusterWeightVector;
        _hypernodeMapping = other._hypernodeMapping;
        _unionFind = UnionFind(_N);
        return *this;
    }

    Graph& operator=(Graph&& other) {
        _N = std::move(other._N);
        _adj_array = std::move(other._adj_array);
        _nodes = std::move(other._nodes);
        _shuffleNodes = std::move(other._shuffleNodes);
        _edges = std::move(other._edges);
        _selfloop = std::move(other._selfloop);
        _weightedDegree = std::move(other._weightedDegree);
        _cluster_id = std::move(other._cluster_id);
        _incidentClusterWeight.assign(_N,IncidentClusterWeight(0,0.0L));
        _total_weight = std::move(other._total_weight);
        _posInIncidentClusterWeightVector = other._posInIncidentClusterWeightVector;
        _hypernodeMapping = std::move(other._hypernodeMapping);   
        _unionFind = UnionFind(_N);
        _posInIncidentClusterWeightVector.clear();
        return *this;
    }
    
    
    
    std::pair<NodeIterator,NodeIterator> nodes() const {
        return std::make_pair(_nodes.begin(),_nodes.end());
    }
    
    std::pair<NodeIterator,NodeIterator> randomNodeOrder() {
        Randomize::instance().shuffleVector(_shuffleNodes, _shuffleNodes.size());
        return std::make_pair(_shuffleNodes.begin(),_shuffleNodes.end());
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
    
    void setHypernodeClusterID(const HypernodeID hn, const ClusterID c_id) {
        ASSERT(_hypernodeMapping[hn] != INVALID_NODE, "Hypernode " << hn << " isn't part of the current graph!");
        _cluster_id[_hypernodeMapping[hn]] = c_id;
    }
    
    void setHyperedgeClusterID(const HyperedgeID he, const ClusterID c_id, const size_t N) {
        ASSERT(_hypernodeMapping[N+he] != INVALID_NODE, "Hyperedge " << he << " isn't part of the current graph!");
        _cluster_id[_hypernodeMapping[N+he]] = c_id;
    }
    
    ClusterID hypernodeClusterID(const HypernodeID hn) const {
        ASSERT(_hypernodeMapping[hn] != INVALID_NODE, "Hypernode " << hn << " isn't part of the current graph!");
        return _cluster_id[_hypernodeMapping[hn]];
    }
    
    
    ClusterID hyperedgeClusterID(const HyperedgeID he, const size_t N) const {
        ASSERT(_hypernodeMapping[N+he] != INVALID_NODE, "Hyperedge " << he << " isn't part of the current graph!");  
        return _cluster_id[_hypernodeMapping[N+he]];
    }
    
    ClusterID clusterID(const NodeID node) const  {
        ASSERT(node < numNodes(), "NodeID " << node << " doesn't exist!");
        return _cluster_id[node];
    }
    
    void setClusterID(const NodeID node, const ClusterID c_id) {
        ASSERT(node < numNodes(), "NodeID " << node << " doesn't exist!");
        _cluster_id[node] = c_id;
    }
    
    void contractHypernodes(const HypernodeID u, const HypernodeID v) {
        //ASSERT(_hypernodeMapping[u] != INVALID_NODE && _hypernodeMapping[v] != INVALID_NODE, "Hypernode " << u << " or " << v << " isn't part of the current graph!");
        NodeID u_id = _hypernodeMapping[u];
        NodeID v_id = _hypernodeMapping[v];
        if(u_id == INVALID_NODE || v_id == INVALID_NODE) return;
        if(u_id != v_id && clusterID(u_id) == clusterID(v_id)) {
            _unionFind.unionSets(u_id,v_id);
        }
    }
    
    Graph contractHyperedgeCluster(NodeID he_begin) {
        _unionFind.reset();
        std::sort(_shuffleNodes.begin(),_shuffleNodes.end());
        std::sort(_shuffleNodes.begin()+he_begin,_shuffleNodes.end(),[&](const NodeID n1, const NodeID n2) {
            return clusterID(n1) < clusterID(n2);
        });
        NodeID cur_node = _shuffleNodes[he_begin];
        ClusterID cur_cid = clusterID(cur_node);
        for(NodeID n = he_begin+1; n < _shuffleNodes.size(); ++n) {
            NodeID node = _shuffleNodes[n];
            ClusterID cid = clusterID(node);
            if(cur_cid == cid) {
                _unionFind.unionSets(cur_node,node);
            }
            else {
                cur_node = node;
                cur_cid = cid;
            }
        }
        return contractGraphWithUnionFind();
    }
    
    Graph contractGraphWithUnionFind() {
        
        std::vector<std::vector<Edge>> adj_list(_N,std::vector<Edge>());
        for(NodeID node : nodes()) {
            NodeID rep = _unionFind.findSet(node);
            for(Edge e : adjacentNodes(node)) {
                Edge new_e;
                new_e.targetNode = _unionFind.findSet(e.targetNode);
                new_e.weight = e.weight;
                adj_list[rep].push_back(new_e);
            }
        }
        
        std::vector<NodeID> mapping(_N,INVALID_NODE);
        NodeID cur_node = 0;
        for(NodeID node = 0; node < _N; ++node) {
            if(_unionFind.findSet(node) == node) {
                mapping[node] = cur_node++;
            }
        }
        
        std::vector<NodeID> hypernodeMapping(_hypernodeMapping.size(),INVALID_NODE);
        for(HypernodeID hn = 0; hn < _hypernodeMapping.size(); ++hn) {
            if(_hypernodeMapping[hn] != INVALID_NODE) {
                hypernodeMapping[hn] = mapping[_hypernodeMapping[hn]];
            }
        }
        
        std::vector<NodeID> adj_array;
        std::vector<Edge> edges;
        std::vector<ClusterID> cluster_id(cur_node,INVALID_CLUSTER);
        std::vector<ClusterID> cluster_id_mapping(_cluster_id.size(),INVALID_CLUSTER);
        ClusterID cur_cid = 0;
        for(NodeID node = 0; node < _N; ++node) {
            if(_unionFind.findSet(node) == node) {
                std::sort(adj_list[node].begin(),adj_list[node].end(),[&](const Edge& e1, const Edge& e2) {
                    return mapping[e1.targetNode] < mapping[e2.targetNode];
                });
                //Adding Sentinel
                Edge sentinel;
                sentinel.targetNode = std::numeric_limits<NodeID>::max();
                sentinel.weight = 0;
                adj_list[node].push_back(sentinel);
                
                ClusterID node_cid = clusterID(node);
                if(cluster_id_mapping[node_cid] == INVALID_CLUSTER) {
                    cluster_id_mapping[node_cid] = cur_cid++;
                }
                cluster_id[mapping[node]] = cluster_id_mapping[node_cid];
                
                adj_array.push_back(edges.size());
                Edge e;
                
                if(adj_list[node][0].targetNode == INVALID_NODE) {
                    continue;
                }
                e.targetNode = mapping[adj_list[node][0].targetNode]; 
                e.weight = adj_list[node][0].weight;;
                for(size_t i = 1; i < adj_list[node].size(); ++i) {
                    NodeID cur_id;
                    
                    if(i < adj_list[node].size()-1) cur_id = mapping[adj_list[node][i].targetNode];
                    else cur_id = adj_list[node][i].targetNode;
                    
                    if(e.targetNode == cur_id) {
                        e.weight += adj_list[node][i].weight;
                    }
                    else {
                        edges.push_back(e);
                        e.targetNode = cur_id;
                        e.weight = adj_list[node][i].weight; 
                    }
                }
            }
        }
        adj_array.push_back(edges.size());
        
        
        Graph graph(adj_array,edges,hypernodeMapping,cluster_id);
        return std::move(graph);
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
    
    void printGraph() {
        
        std::cout << "Number Nodes: " << numNodes() << std::endl;
        std::cout << "Number Edges: " << numEdges() << std::endl;
        
        for(NodeID n : nodes()) {
            std::cout << "Node ID: " << n << "(Comm.: "<<clusterID(n)<<"), Adj. List: ";
            for(Edge e : adjacentNodes(n)) {
                std::cout << "(" << e.targetNode << ",w=" << e.weight<<") ";
            }
            std::cout << "\n";
        }
        
    }
    
protected:
    
    NodeID _N;
    std::vector<NodeID> _adj_array;
    std::vector<NodeID> _nodes;
    std::vector<NodeID> _shuffleNodes;
    std::vector<Edge> _edges;
    std::vector<EdgeWeight> _selfloop;
    std::vector<EdgeWeight> _weightedDegree;
    std::vector<ClusterID> _cluster_id;
    std::vector<IncidentClusterWeight> _incidentClusterWeight;
    EdgeWeight _total_weight;
    SparseMap<ClusterID,size_t> _posInIncidentClusterWeightVector;
    std::vector<NodeID> _hypernodeMapping;
    UnionFind _unionFind;
    
private:
    //FRIEND_TEST(AGraph,DeterminesIncidentClusterWeightsOfAClusterCorrect);
    
    
    /**
     * Creates an iterator to all incident Clusters of ClusterID cid. Iterator points to an 
     * IncidentClusterWeight-Struct which contains the incident ClusterID clusterID and the sum of
     * the weights of all incident edges from cid to clusterID.
     * 
     * @param cid ClusterID, which incident clusters should be evaluated
     * @return Iterator to all incident clusters of ClusterID cid
     */
    std::pair<IncidentClusterWeightIterator,IncidentClusterWeightIterator> incidentClusterWeightOfCluster(std::pair<NodeIterator,NodeIterator>& cluster_range);
    
    
    void constructBipartiteGraph(const Hypergraph& hg, const bool use_uniform_edge_weight) {
        NodeID sum_edges = 0;
        
        size_t N = static_cast<size_t>(hg.initialNumNodes());
        
        NodeID cur_node_id = 0;
        
        //Construct adj. array for all hypernode. Amount of edges is equal to the degree of the corresponding hypernode.
        for(HypernodeID hn : hg.nodes()) {
            _hypernodeMapping[hn] = cur_node_id;
            _adj_array[cur_node_id++] = sum_edges;
            sum_edges += hg.nodeDegree(hn);
        }
        
        //Construct adj. array for all hyperedges. Amount of edges is equal to the size of the corresponding hyperedge.
        for(HyperedgeID he : hg.edges()) {
            _hypernodeMapping[N+he] = cur_node_id;
            _adj_array[cur_node_id++] = sum_edges;
            sum_edges += hg.edgeSize(he);
        }
        
        _adj_array[_N] = sum_edges;
        _edges.resize(sum_edges);
        
        for(HypernodeID hn : hg.nodes()) {
            
            size_t pos = 0;
            for(HyperedgeID he : hg.incidentEdges(hn)) {
                Edge e;
                e.targetNode = _hypernodeMapping[N + he];
                if(!use_uniform_edge_weight) {
                    e.weight = static_cast<EdgeWeight>(hg.edgeWeight(he))/
                               static_cast<EdgeWeight>(hg.edgeSize(he));
                }
                else {
                    e.weight = static_cast<EdgeWeight>(hg.edgeWeight(he));
                }
                _total_weight += e.weight;
                _weightedDegree[_hypernodeMapping[hn]] += e.weight;
                _edges[_adj_array[_hypernodeMapping[hn]] + pos++] = e;
            }
        }
        
        for(HyperedgeID he : hg.edges()) {
           size_t pos = 0;
           for(HypernodeID hn : hg.pins(he)) {
                Edge e;
                e.targetNode = _hypernodeMapping[hn];
                if(!use_uniform_edge_weight) {
                    e.weight = static_cast<EdgeWeight>(hg.edgeWeight(he))/
                    static_cast<EdgeWeight>(hg.edgeSize(he));
                }
                else {
                    e.weight = static_cast<EdgeWeight>(hg.edgeWeight(he));
                }
                _total_weight += e.weight;
                _weightedDegree[_hypernodeMapping[N + he]] += e.weight;
                _edges[_adj_array[_hypernodeMapping[N + he]]+pos++] = e;
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

std::pair<IncidentClusterWeightIterator,IncidentClusterWeightIterator> Graph::incidentClusterWeightOfCluster(std::pair<NodeIterator,NodeIterator>& cluster_range) {
    
    _posInIncidentClusterWeightVector.clear();
    size_t idx = 0;

    for(NodeID node : cluster_range) {
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
    
    
    std::sort(_shuffleNodes.begin(),_shuffleNodes.end());
    std::sort(_shuffleNodes.begin(),_shuffleNodes.end(),[&](const NodeID& n1, const NodeID& n2) {
        return _cluster_id[n1] < _cluster_id[n2]  || (_cluster_id[n1] == _cluster_id[n2] && n1 < n2);
    });
    
    //Add Sentinels
    _shuffleNodes.push_back(_cluster_id.size());
    _cluster_id.push_back(new_cid);
    
    std::vector<NodeID> new_adj_array(new_cid+1,0);
    std::vector<Edge> new_edges;
    size_t start_idx = 0;
    for(size_t i = 0; i < _N+1; ++i) {
        if(_cluster_id[_shuffleNodes[start_idx]] != _cluster_id[_shuffleNodes[i]]) {
            ClusterID cid = _cluster_id[_shuffleNodes[start_idx]];
            new_adj_array[cid] = new_edges.size();
            std::pair<NodeIterator,NodeIterator> cluster_range = std::make_pair(_shuffleNodes.begin()+start_idx,
                                                                                _shuffleNodes.begin()+i);
            for(auto incidentClusterWeight : incidentClusterWeightOfCluster(cluster_range)) {
                Edge e;
                e.targetNode = static_cast<NodeID>(incidentClusterWeight.clusterID);
                e.weight = incidentClusterWeight.weight;
                new_edges.push_back(e);
            }
            start_idx = i;
        }
    }
    
    //Remove Sentinels
    _shuffleNodes.pop_back(); _cluster_id.pop_back();
    
    new_adj_array[new_cid] = new_edges.size();
    Graph graph(new_adj_array,new_edges);
    
    return std::make_pair(std::move(graph),node2contractedNode);
}
        
}  // namespace ds
}  // namespace kahypar

