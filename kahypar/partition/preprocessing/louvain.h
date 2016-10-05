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
using ds::Edge;

template<class QualityMeasure = Mandatory>
class Louvain {
    
public:
    
    Louvain(const Hypergraph& hypergraph) : graph(hypergraph) { }
    
    
    void louvain() {
        bool improvement = false;
        
        std::vector<std::vector<NodeID>> mapping_stack;
        std::vector<Graph> graph_stack;
        
        graph_stack.push_back(graph);
        int cur_idx = 0;
        
        do {
            QualityMeasure quality(graph_stack[cur_idx]);
            
            improvement = louvain_pass(graph_stack[cur_idx],quality);
            
            if(improvement) {
                auto contraction = graph_stack[cur_idx++].contractCluster();
                graph_stack.push_back(std::move(contraction.first));
                mapping_stack.push_back(std::move(contraction.second));
            }
            
        } while(improvement);
        
        while(!mapping_stack.empty()) {
            assignClusterToNextLevelFinerGraph(graph_stack[cur_idx-1],graph_stack[cur_idx],mapping_stack[cur_idx-1]);
            graph_stack.pop_back();
            mapping_stack.pop_back();
        }
        
        for(NodeID node : graph.nodes()) {
            graph.setClusterID(node,graph_stack[0].clusterID(node));    
        }    
    }

private:
    FRIEND_TEST(ALouvainAlgorithm,DoOneLouvainPass);
    FRIEND_TEST(ALouvainAlgorithm,AssingsMappingToNextLevelFinerGraph);
    
    void assignClusterToNextLevelFinerGraph(Graph& fineGraph, Graph& coarseGraph, std::vector<NodeID>& mapping) {
        for(NodeID node : fineGraph.nodes()) {
            fineGraph.setClusterID(node,coarseGraph.clusterID(mapping[node]));
        }
    }
    
    bool louvain_pass(Graph& g, QualityMeasure& quality) {
        bool improvement = false;
        size_t node_moves = 0;
        EdgeWeight cur_quality = quality.quality();
        EdgeWeight new_quality = cur_quality;
        
        //TODO(heuer): Think about shuffling nodes before louvain pass

        do {
            node_moves = 0;
            
            for(NodeID node : g.nodes()) {
                ClusterID cur_cid = g.clusterID(node);
                ClusterID best_cid = cur_cid;
                EdgeWeight best_incident_cluster_weight = 0.0L;
                EdgeWeight best_gain = 0.0L;
                
                for(Edge e : g.adjacentNodes(node)) {
                    if(g.clusterID(e.targetNode) == cur_cid) {
                        best_incident_cluster_weight += e.weight;
                    }
                }
                
                quality.remove(node,best_incident_cluster_weight);
                
                for(auto cluster : g.incidentClusterWeightOfNode(node)) {
                    ClusterID cid = cluster.clusterID;
                    EdgeWeight weight = cluster.weight;
                    EdgeWeight gain = quality.gain(node,cid,weight);
                    if(gain > best_gain) {
                        best_gain = gain;
                        best_incident_cluster_weight = weight;
                        best_cid = cid;
                    }
                }
                
                quality.insert(node,best_cid,best_incident_cluster_weight);
                
                if(best_cid != cur_cid) {
                    node_moves++;
                }
                
            }
            
        } while(node_moves > 0);
        
        new_quality = quality.quality();
        improvement = (new_quality - cur_quality) > 0.0L;

        return improvement;
    }
    
    Graph graph;
};

}  // namespace kahypar

