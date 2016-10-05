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

#define EPS 1e-5

template<class QualityMeasure = Mandatory>
class Louvain {
    
public:
    
    Louvain(const Hypergraph& hypergraph) : graph(hypergraph), _num_hypernodes(hypergraph.initialNumNodes()) { }
    
    
    void louvain() {
        bool improvement = false;
        EdgeWeight old_quality = -1.0L;
        EdgeWeight cur_quality = -1.0L;
        
        std::vector<std::vector<NodeID>> mapping_stack;
        std::vector<Graph> graph_stack;
        
        graph_stack.push_back(graph);
        int cur_idx = 0;
        
        do {
            QualityMeasure quality(graph_stack[cur_idx]);
            
            //Checks if quality of the coarse graph is equal with the quality of next level finer graph
            ASSERT([&]() {
                if(cur_idx == 0) return true;
                if(std::abs(cur_quality-quality.quality()) > EPS) {
                    LOGVAR(cur_quality);
                    LOGVAR(quality.quality());
                    return false;
                }
                return true;
            }(),"Quality of the contracted graph is not equal with its corresponding uncontracted graph!");
            
            old_quality = cur_quality;
            cur_quality = louvain_pass(graph_stack[cur_idx],quality);
            improvement = cur_quality - old_quality > EPS;
            
            LOG("Old Quality: " << old_quality << " - New Quality: " << cur_quality);
            
            if(improvement) {
                ASSERT(cur_quality-old_quality >= EPS, "Quality should be improved during louvain pass!");
                cur_quality = quality.quality();
                auto contraction = graph_stack[cur_idx++].contractCluster();
                graph_stack.push_back(std::move(contraction.first));
                mapping_stack.push_back(std::move(contraction.second));
            }
            
        } while(improvement);
        
        ASSERT((mapping_stack.size() + 1) == graph_stack.size(), "Unequality between graph and mapping stack!");
        while(!mapping_stack.empty()) {
            assignClusterToNextLevelFinerGraph(graph_stack[cur_idx-1],graph_stack[cur_idx],mapping_stack[cur_idx-1]);
            graph_stack.pop_back();
            mapping_stack.pop_back();
            cur_idx--;
        }
        
        for(NodeID node : graph.nodes()) { 
            graph.setClusterID(node,graph_stack[0].clusterID(node));    
        }    
        
    }
    
    std::vector<ClusterID> getClusterIDsForAllHypernodes() {
        std::vector<ClusterID> clusterIDs(_num_hypernodes,0);
        for(NodeID node : graph.nodes()) {
            clusterIDs[node] = graph.clusterID(node);
        }
        return clusterIDs;
    }

private:
    FRIEND_TEST(ALouvainAlgorithm,DoesOneLouvainPass);
    FRIEND_TEST(ALouvainAlgorithm,AssingsMappingToNextLevelFinerGraph);
    
    void assignClusterToNextLevelFinerGraph(Graph& fineGraph, Graph& coarseGraph, std::vector<NodeID>& mapping) {
        for(NodeID node : fineGraph.nodes()) {
            fineGraph.setClusterID(node,coarseGraph.clusterID(mapping[node]));
        }
        
        //Check if cluster ids from coarse graph are succesfully mapped to the nodes of 
        //next level finer graph
        ASSERT([&]() {
            for(NodeID node : fineGraph.nodes()) {
                if(fineGraph.clusterID(node) != coarseGraph.clusterID(mapping[node])) {
                    LOGVAR(fineGraph.clusterID(node));
                    LOGVAR(coarseGraph.clusterID(mapping[node]));
                    return false;
                }
            }
            return true;
        }(),"Mapping from coarse to fine graph failed!");
    }
    
    EdgeWeight louvain_pass(Graph& g, QualityMeasure& quality) {
        size_t node_moves = 0;
        
        //TODO(heuer): Think about shuffling nodes before louvain pass

        do {
            
            node_moves = 0;
            
            for(NodeID node : g.nodes()) {
                ClusterID cur_cid = g.clusterID(node);
                EdgeWeight cur_incident_cluster_weight = 0.0L;
                ClusterID best_cid = cur_cid;
                EdgeWeight best_incident_cluster_weight = 0.0L;
                EdgeWeight best_gain = 0.0L;
                
                for(Edge e : g.adjacentNodes(node)) {
                    if(g.clusterID(e.targetNode) == cur_cid && e.targetNode != node) {
                        cur_incident_cluster_weight += e.weight;
                    }
                }
                best_incident_cluster_weight = cur_incident_cluster_weight;
                
                quality.remove(node,cur_incident_cluster_weight);
                
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
                    ASSERT([&]() {
                        quality.remove(node,best_incident_cluster_weight); //Remove node from best cluster...
                        quality.insert(node,cur_cid,cur_incident_cluster_weight); // ... and insert in his old cluster.
                        EdgeWeight quality_before = quality.quality();
                        quality.remove(node,cur_incident_cluster_weight); //Remove node again from his old cluster ...
                        quality.insert(node,best_cid,best_incident_cluster_weight); //... and insert it in cluster with best gain.
                        EdgeWeight quality_after = quality.quality();
                        if(quality_after - quality_before < -EPS) {
                            LOGVAR(quality_before);
                            LOGVAR(quality_after);
                            return false;
                        }
                        return true;
                    }(),"Move did not increase the quality!");
                    node_moves++;
                }
                
            }
            
        } while(node_moves > 0);
        

        return quality.quality();
    }
    
    Graph graph;
    HypernodeID _num_hypernodes;
};

}  // namespace kahypar

