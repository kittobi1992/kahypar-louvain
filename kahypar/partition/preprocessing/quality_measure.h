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

namespace kahypar {

    
class QualityMeasure {

friend class Modularity;
    
public:
    
    QualityMeasure(Graph& graph) : graph(graph) { }
    
    virtual ~QualityMeasure() { }
    
    virtual void remove(NodeID node, EdgeWeight incidentCommWeight)=0;
    virtual void insert(NodeID node, ClusterID new_cid, EdgeWeight incidentCommWeight)=0;
    virtual EdgeWeight gain(NodeID node, ClusterID cid, EdgeWeight incidentCommWeight)=0;
    virtual EdgeWeight quality()=0;

private:
    Graph& graph;
    
};

class Modularity : public QualityMeasure {

public:
    Modularity(Graph& graph) : QualityMeasure(graph), in(graph.numNodes()), tot(graph.numNodes()) { 
        for(NodeID node : graph.nodes()) {
            in[node] = graph.selfloopWeight(node);
            tot[node] = graph.weightedDegree(node);
        }
    }
    
    ~Modularity() {
        in.clear();
        tot.clear();
    }
    
    inline void remove(NodeID node, EdgeWeight incidentCommWeight) {
        ClusterID cid = graph.clusterID(node);
        
        in[cid] -= 2.0L*incidentCommWeight + graph.selfloopWeight(node);
        tot[cid] -= graph.weightedDegree(node);
        
        graph.setClusterID(node,-1);
    }
    
    inline void insert(NodeID node, ClusterID new_cid, EdgeWeight incidentCommWeight) {
        
        in[new_cid] += 2.0L*incidentCommWeight + graph.selfloopWeight(node);
        tot[new_cid] += graph.weightedDegree(node);
        
        graph.setClusterID(node,new_cid);
    }
    
    inline EdgeWeight gain(NodeID node, ClusterID cid, EdgeWeight incidentCommWeight) {
        
        EdgeWeight totc = tot[cid];
        EdgeWeight m2 = graph.totalWeight();
        EdgeWeight w_degree = graph.weightedDegree(node);
        
        return (incidentCommWeight - totc*w_degree/m2);
    }
    
    
    EdgeWeight quality() {
        EdgeWeight q = 0.0L;
        EdgeWeight m2 = graph.totalWeight();
        
        for(NodeID node : graph.nodes()) {
            if(tot[node] > 0.0L) {
                q += in[node] - (tot[node]*tot[node])/m2;
            }
        }
        
        q /= m2;
        
        return q;
    }
    
private:
    FRIEND_TEST(AModularityMeasure,IsCorrectInitialized);    FRIEND_TEST(AModularityMeasure,RemoveNodeFromCommunity);
    FRIEND_TEST(AModularityMeasure,InsertNodeInCommunity);
    FRIEND_TEST(AModularityMeasure,RemoveNodeFromCommunityWithMoreThanOneNode);
    
    using QualityMeasure::graph;
    std::vector<EdgeWeight> in;
    std::vector<EdgeWeight> tot;
    
    
};

}  // namespace kahypar

