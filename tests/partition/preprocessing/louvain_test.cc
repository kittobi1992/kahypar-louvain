/***************************************************************************
 *  Copyright (C) 2015 Tobias Heuer <tobias.heuer@gmx.net>
 **************************************************************************/

#include <vector>
#include <set>

#include "gmock/gmock.h"

#include "kahypar/definitions.h"
#include "kahypar/macros.h"
#include "kahypar/partition/preprocessing/louvain.h"
#include "kahypar/partition/preprocessing/quality_measure.h"

using::testing::Eq;
using::testing::Test;


namespace kahypar {
    
#define INVALID -1
#define EPS 1e-5
   
class ALouvainAlgorithm : public Test {
public:
    ALouvainAlgorithm() : louvain(nullptr),
                        hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9, 12 },
                                   HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 }) { 
        louvain = std::make_shared<Louvain<Modularity>>(hypergraph);
    }
                   
    std::shared_ptr<Louvain<Modularity>> louvain;
    Hypergraph hypergraph;
};

class AModularityMeasure : public Test {
public:
    AModularityMeasure() : modularity(nullptr),
                           hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9, 12 },
                                 HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 }),
                           graph(hypergraph) { 
        modularity = std::make_shared<Modularity>(graph);
    }
               
    std::shared_ptr<Modularity> modularity;
    Hypergraph hypergraph;    
    Graph graph;
};

TEST_F(AModularityMeasure,IsCorrectInitialized) {
    std::vector<EdgeWeight> expected_in = {0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L};
    std::vector<EdgeWeight> expected_tot = {0.75L,0.25L,0.5L+1.0L/3.0L,0.25L+1.0L/3.0L,0.25L+1.0L/3.0L,1.0L/3.0L,2.0L/3.0L,1.0L,1.0L,1.0L,1.0L};
    for(NodeID node : graph.nodes()) {
        ASSERT_LE(std::abs(expected_in[node]-modularity->in[node]),EPS);
        ASSERT_LE(std::abs(expected_tot[node]-modularity->tot[node]),EPS);
    }
}

TEST_F(AModularityMeasure,RemoveNodeFromCommunity) {
    for(NodeID node : graph.nodes()) {
        modularity->remove(node,0.0L);
        ASSERT_EQ(0.0L,modularity->in[node]);
        ASSERT_EQ(0.0L,modularity->tot[node]);
        ASSERT_EQ(graph.clusterID(node),-1);
    }
}

TEST_F(AModularityMeasure,RemoveNodeFromCommunityWithMoreThanOneNode) {
    modularity->remove(1,0.0L);
    modularity->insert(1,8,0.25L);
    modularity->remove(3,0.0L);
    modularity->insert(3,8,0.25L);
    
    modularity->remove(8,0.5);
    
    ASSERT_EQ(0.0L,modularity->in[8]);
    ASSERT_LE(std::abs(0.5L+1.0L/3.0L-modularity->tot[8]),EPS);
    ASSERT_EQ(graph.clusterID(8),-1);
}

TEST_F(AModularityMeasure,InsertNodeInCommunity) {
    modularity->remove(1,0.0L);
    modularity->insert(1,8,0.25L);
    ASSERT_EQ(0.0L,modularity->in[1]);
    ASSERT_EQ(0.0L,modularity->tot[1]);
    ASSERT_EQ(0.5L,modularity->in[8]);
    ASSERT_EQ(1.25L,modularity->tot[8]);
    ASSERT_EQ(graph.clusterID(1),8);
    
    modularity->remove(3,0.0L);
    modularity->insert(3,8,0.25L);
    ASSERT_EQ(0.0L,modularity->in[3]);
    ASSERT_EQ(0.0L,modularity->tot[3]);
    ASSERT_EQ(1.0L,modularity->in[8]);
    ASSERT_LE(std::abs(1.25L+(0.25L+1.0L/3.0L)-modularity->tot[8]),EPS);
    ASSERT_EQ(graph.clusterID(3),8);
}

TEST_F(AModularityMeasure,CalculatesCorrectGainValuesForIsolatedNode) {
    modularity->remove(1,0.0L);
    modularity->insert(1,8,0.25L);
    modularity->remove(3,0.0L);
    modularity->insert(3,8,0.25L);
    
    modularity->remove(4,0.0);
    
    for(auto incidentClusterWeight : graph.incidentClusterWeightOfNode(4)) {
        ClusterID cid = incidentClusterWeight.clusterID;
        EdgeWeight w = incidentClusterWeight.weight;
        EdgeWeight gain = modularity->gain(4,cid,w);
        EdgeWeight expected_gain = std::numeric_limits<EdgeWeight>::max()/2.0L;
        if(cid == 8) expected_gain = 0.116319444L;
        else if(cid == 9) expected_gain = 0.260416667L;
        ASSERT_LE(std::abs(expected_gain-gain),EPS);
    }
}

TEST_F(ALouvainAlgorithm,DoOneLouvainPass) { 
    Graph graph(hypergraph);
    Modularity modularity(graph);
    EdgeWeight quality_before = modularity.quality();
    louvain->louvain_pass(graph,modularity);
    EdgeWeight quality_after = modularity.quality();
    ASSERT_LE(quality_before,quality_after);
}

TEST_F(ALouvainAlgorithm,AssingsMappingToNextLevelFinerGraph) {
    Graph graph(hypergraph);
    Modularity modularity(graph);
    louvain->louvain_pass(graph,modularity);
    auto contraction = graph.contractCluster();
    Graph coarseGraph = std::move(contraction.first);
    std::vector<NodeID> mapping = std::move(contraction.second);
    louvain->assignClusterToNextLevelFinerGraph(graph,coarseGraph,mapping);
    ASSERT_EQ(0,graph.clusterID(0));
    ASSERT_EQ(1,graph.clusterID(1));
    ASSERT_EQ(0,graph.clusterID(2));
    ASSERT_EQ(1,graph.clusterID(3));
    ASSERT_EQ(1,graph.clusterID(4));
    ASSERT_EQ(2,graph.clusterID(5));
    ASSERT_EQ(2,graph.clusterID(6));
    ASSERT_EQ(0,graph.clusterID(7));
    ASSERT_EQ(1,graph.clusterID(8));
    ASSERT_EQ(1,graph.clusterID(9));
    ASSERT_EQ(2,graph.clusterID(10));
}

TEST_F(ALouvainAlgorithm,DoLouvainAlgorithm) { 
    louvain->louvain();
}

} //namespace kahypar