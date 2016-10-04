/***************************************************************************
 *  Copyright (C) 2015 Tobias Heuer <tobias.heuer@gmx.net>
 **************************************************************************/

#include <vector>
#include <set>

#include "gmock/gmock.h"

#include "kahypar/definitions.h"
#include "kahypar/datastructure/hgr_to_graph.h"

using::testing::Eq;
using::testing::Test;


namespace kahypar {

#define INVALID -1
#define EPS 1e-5
    
using NodeID = HypernodeID;
using EdgeWeight = long double;
using ClusterID = PartitionID;
using ds::Edge;
using ds::IncidentClusterWeight;
    
class AHgrToGraphTest : public Test {
public:
    AHgrToGraphTest() : graph(nullptr),
                        hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9, 12 },
                                   HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 }) { 
        graph = std::make_shared<ds::HgrToGraph>(hypergraph);
    }
                   
    std::shared_ptr<ds::HgrToGraph> graph;
    Hypergraph hypergraph;
};

bool bicoloring(NodeID cur, int cur_col, std::vector<int>& col, std::shared_ptr<ds::HgrToGraph>& graph) {
    col[cur] = cur_col;
    for(Edge e : graph->adjacentNodes(cur)) {
        NodeID id = e.targetNode;
        if(col[id] == INVALID) return bicoloring(id,1-cur_col,col,graph);
        else if(col[id] == col[cur]) return false;
    }
    return true;
}


TEST_F(AHgrToGraphTest, ChecksIfGraphIsBipartite) {
    bool isBipartite = true;
    std::vector<int> col(graph->numNodes(),INVALID);
    for(NodeID id : graph->nodes()) {
        if(col[id] == INVALID) isBipartite && bicoloring(id,0,col,graph);
    }
    ASSERT_TRUE(isBipartite);
}


TEST_F(AHgrToGraphTest, ChecksCorrectConstruction) {
    HypernodeID numNodes = hypergraph.initialNumNodes();
    for(HypernodeID hn : hypergraph.nodes()) {
        std::set<NodeID> edges;
        for(HyperedgeID he : hypergraph.incidentEdges(hn)) {
            edges.insert(numNodes+he);
        }
        ASSERT_EQ(edges.size(),graph->degree(hn));
        for(Edge e : graph->adjacentNodes(hn)) {
            ASSERT_FALSE(edges.find(e.targetNode) == edges.end());
        }
    }
    
    for(HyperedgeID he : hypergraph.edges()) {
        std::set<NodeID> edges;
        for(HypernodeID pin : hypergraph.pins(he)) {
            edges.insert(pin);
        }
        ASSERT_EQ(edges.size(),graph->degree(he+numNodes));
        for(Edge e : graph->adjacentNodes(he+numNodes)) {
            ASSERT_FALSE(edges.find(e.targetNode) == edges.end());
            ASSERT_EQ(e.weight,static_cast<EdgeWeight>(hypergraph.edgeWeight(he))/
                               static_cast<EdgeWeight>(hypergraph.edgeSize(he)));
        }
    }
}

TEST_F(AHgrToGraphTest, ChecksIfClusterIDsAreInitializedCorrectly) {
    for(NodeID id : graph->nodes()) {
        ASSERT_EQ(id,graph->clusterID(id));
    }
}

TEST_F(AHgrToGraphTest, ChecksIfIncidentClusterWeightsIsDeterminedCorrect) {
    NodeID node = 8;
    std::vector<EdgeWeight> cluster_weight = {0.25L,0.25L,0.0L,0.25L,0.25L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    std::vector<bool> incident_cluster = {true,true,false,true,true,false,false,false,true,false,false};
    for(IncidentClusterWeight incidentClusterWeight : graph->incidentClusterWeight(node)) {
        ClusterID c_id = incidentClusterWeight.clusterID;
        EdgeWeight weight = incidentClusterWeight.weight;
        ASSERT_TRUE(incident_cluster[c_id]);
        ASSERT_EQ(cluster_weight[c_id],weight);
    }
}

TEST_F(AHgrToGraphTest, ChecksIfIncidentClusterWeightsIsDeterminedCorrect2) {
    NodeID node = 8;
    graph->setClusterID(1,0);
    graph->setClusterID(3,0);
    std::vector<EdgeWeight> cluster_weight = {0.75L,0.0L,0.0L,0.0L,0.25L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    std::vector<bool> incident_cluster = {true,false,false,false,true,false,false,false,true,false,false};
    for(IncidentClusterWeight incidentClusterWeight : graph->incidentClusterWeight(node)) {
        ClusterID c_id = incidentClusterWeight.clusterID;
        EdgeWeight weight = incidentClusterWeight.weight;
        ASSERT_TRUE(incident_cluster[c_id]);
        ASSERT_EQ(cluster_weight[c_id],weight);
    }
}

TEST_F(AHgrToGraphTest, ChecksIfIncidentNeighborClustersAreDeterminedCorrect) {
    graph->setClusterID(2,0);
    graph->setClusterID(7,0);
    graph->setClusterID(1,1);
    graph->setClusterID(3,1);
    graph->setClusterID(4,1);
    graph->setClusterID(8,1);
    graph->setClusterID(9,1);
    graph->setClusterID(5,2);
    graph->setClusterID(6,2);
    graph->setClusterID(10,2);
    std::vector<EdgeWeight> cluster_weight = {1.0L,0.25L,1.0L/3.0L};
    std::vector<bool> incident_cluster = {true,true,true};
    for(auto incidentClusterWeight : graph->incidentNeighborClusters(0)) {
        ClusterID c_id = incidentClusterWeight.clusterID;
        EdgeWeight weight = incidentClusterWeight.weight; 
        ASSERT_TRUE(incident_cluster[c_id]);
        ASSERT_LE(std::abs(cluster_weight[c_id]-weight),EPS);
    }
    cluster_weight = {1.0L/4.0L,0.75L+2.0L/3.0L,1.0L/3.0L};
    incident_cluster = {true,true,true};
    for(auto incidentClusterWeight : graph->incidentNeighborClusters(1)) {
        ClusterID c_id = incidentClusterWeight.clusterID;
        EdgeWeight weight = incidentClusterWeight.weight; 
        ASSERT_TRUE(incident_cluster[c_id]);
        ASSERT_LE(std::abs(cluster_weight[c_id]-weight),EPS);
    }
    cluster_weight = {1.0L/3.0L,1.0L/3.0L,2.0L/3.0L};
    incident_cluster = {true,true,true};
    for(auto incidentClusterWeight : graph->incidentNeighborClusters(2)) {
        ClusterID c_id = incidentClusterWeight.clusterID;
        EdgeWeight weight = incidentClusterWeight.weight; 
        ASSERT_TRUE(incident_cluster[c_id]);
        ASSERT_LE(std::abs(cluster_weight[c_id]-weight),EPS);
    }
}

TEST_F(AHgrToGraphTest, ChecksIfMappingOfContractedGraphIsCorrect) {
    graph->setClusterID(2,0);
    graph->setClusterID(7,0);
    graph->setClusterID(1,3);
    graph->setClusterID(4,3);
    graph->setClusterID(8,3);
    graph->setClusterID(9,3);
    graph->setClusterID(5,6);
    graph->setClusterID(10,6);
    auto contractedGraph = graph->contractCluster();
    ds::HgrToGraph graph = std::move(contractedGraph.first);
    std::vector<NodeID> mappingToOriginalGraph = contractedGraph.second;
    std::vector<NodeID> correctMappingToOriginalGraph = {0,1,0,1,1,2,2,0,1,1,2};
    for(size_t i = 0; i < mappingToOriginalGraph.size(); ++i) {
        ASSERT_EQ(mappingToOriginalGraph[i],correctMappingToOriginalGraph[i]);
    }
}

TEST_F(AHgrToGraphTest, ChecksIfContractedGraphIsCorrect) {
    graph->setClusterID(2,0);
    graph->setClusterID(7,0);
    graph->setClusterID(1,3);
    graph->setClusterID(4,3);
    graph->setClusterID(8,3);
    graph->setClusterID(9,3);
    graph->setClusterID(5,6);
    graph->setClusterID(10,6);
    auto contractedGraph = graph->contractCluster();
    ds::HgrToGraph graph = std::move(contractedGraph.first);
    std::vector<NodeID> mappingToOriginalGraph = contractedGraph.second;
    
    ASSERT_EQ(3,graph.numNodes());
    std::vector<EdgeWeight> edge_weight = {1.0L,0.25L,1.0L/3.0L};
    std::vector<bool> incident_nodes = {true,true,true};
    for(Edge e : graph.adjacentNodes(0)) {
        NodeID n_id = e.targetNode;
        EdgeWeight weight = e.weight; 
        ASSERT_TRUE(incident_nodes[n_id]);
        ASSERT_LE(std::abs(edge_weight[n_id]-weight),EPS);           
    }
    
    edge_weight = {1.0L/4.0L,0.75L+2.0L/3.0L,1.0L/3.0L};
    incident_nodes = {true,true,true};
    for(Edge e : graph.adjacentNodes(1)) {
        NodeID n_id = e.targetNode;
        EdgeWeight weight = e.weight; 
        ASSERT_TRUE(incident_nodes[n_id]);
        ASSERT_LE(std::abs(edge_weight[n_id]-weight),EPS);           
    }
    
    edge_weight = {1.0L/3.0L,1.0L/3.0L,2.0L/3.0L};
    incident_nodes = {true,true,true};
    for(Edge e : graph.adjacentNodes(2)) {
        NodeID n_id = e.targetNode;
        EdgeWeight weight = e.weight; 
        ASSERT_TRUE(incident_nodes[n_id]);
        ASSERT_LE(std::abs(edge_weight[n_id]-weight),EPS);           
    }
    
}




}