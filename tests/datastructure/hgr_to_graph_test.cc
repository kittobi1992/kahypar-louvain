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
    
using NodeID = HypernodeID;
using EdgeWeight = long double;
using ClusterID = PartitionID;
    
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
    for(NodeID id : graph->edges(cur)) {
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
        for(NodeID id : graph->edges(hn)) {
            ASSERT_FALSE(edges.find(id) == edges.end());
        }
    }
    
    for(HyperedgeID he : hypergraph.edges()) {
        std::set<NodeID> edges;
        for(HypernodeID pin : hypergraph.pins(he)) {
            edges.insert(pin);
        }
        ASSERT_EQ(edges.size(),graph->degree(he+numNodes));
        for(NodeID id : graph->edges(he+numNodes)) {
            ASSERT_FALSE(edges.find(id) == edges.end());
        }
    }
}

TEST_F(AHgrToGraphTest, ChecksIfClusterIDsAreInitializedCorrectly) {
    for(NodeID id : graph->nodes()) {
        ASSERT_EQ(id,graph->clusterID(id));
    }
}

TEST_F(AHgrToGraphTest, ChecksIfIncidentClustersAndWeightsAreDeterminedCorrect) {
    NodeID node = 8;
    std::vector<EdgeWeight> cluster_weight = {0.25L,0.25L,0.0L,0.25L,0.25L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    std::vector<bool> incident_cluster = {true,true,false,true,true,false,false,false,true,false,false};
    auto it = graph->incidentClusterWeight(node);
    for(auto cur = it.first; cur != it.second; ++cur) {
        ClusterID c_id = cur->first;
        EdgeWeight weight = cur->second;
        ASSERT_TRUE(incident_cluster[c_id]);
        ASSERT_EQ(cluster_weight[c_id],weight);
    }
}

TEST_F(AHgrToGraphTest, ChecksIfIncidentClustersAndWeightsAreDeterminedCorrect2) {
    NodeID node = 8;
    graph->setClusterID(1,0);
    graph->setClusterID(3,0);
    std::vector<EdgeWeight> cluster_weight = {0.75L,0.0L,0.0L,0.0L,0.25L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    std::vector<bool> incident_cluster = {true,false,false,false,true,false,false,false,true,false,false};
    auto it = graph->incidentClusterWeight(node);
    for(auto cur = it.first; cur != it.second; ++cur) {
        ClusterID c_id = cur->first;
        EdgeWeight weight = cur->second;
        ASSERT_TRUE(incident_cluster[c_id]);
        ASSERT_EQ(cluster_weight[c_id],weight);
    }
}



}