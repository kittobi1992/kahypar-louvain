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
#include "kahypar/partition/configuration.h"

namespace kahypar {

using ds::Graph;
using ds::Edge;

#define EPS 1e-5

template<class QualityMeasure = Mandatory>
class Louvain {
    
public:
    
    Louvain(const Hypergraph& hypergraph, const Configuration& config) : _graph(hypergraph,config), 
                                                                         _config(config), _first_louvain_call(true) { }
                                                                         
    Louvain(Graph& graph, const Configuration& config) : _graph(graph), _config(config), _first_louvain_call(true) { }    
    
    Louvain(Graph&& graph, const Configuration& config) : _graph(std::move(graph)), _config(config), _first_louvain_call(true) { } 
    
    
    
    EdgeWeight louvain() {
        bool improvement = false;
        size_t iteration = 0;
        EdgeWeight old_quality = -1.0L;
        EdgeWeight cur_quality = -1.0L;
        
        size_t max_iterations = std::numeric_limits<size_t>::max();
        
        if(_config.preprocessing.use_multilevel_louvain) {
            if(_first_louvain_call) {
                max_iterations = 1;
            }
            else {
                if(_config.preprocessing.louvain_contract_graph_like_hg) {
                    max_iterations = 2;
                    _graph = _graph.contractGraphWithUnionFind();
                }
                else {
                    max_iterations = 1;
                    _graph = _graph.contractCluster().first;
                }
            }
        }
        
        std::vector<std::vector<NodeID>> mapping_stack;
        std::vector<Graph> graph_stack;
        
        graph_stack.push_back(_graph);
        int cur_idx = 0;

        do {
            
            LOG("Graph Number Nodes: " << graph_stack[cur_idx].numNodes());
            LOG("Graph Number Edges: " << graph_stack[cur_idx].numEdges());
            QualityMeasure quality(graph_stack[cur_idx]);
            if(iteration == 0) {
                cur_quality = quality.quality();
            }
             
            LOG("######## Starting Louvain-Pass #" << ++iteration << " ########");
            
            //Checks if quality of the coarse graph is equal with the quality of next level finer graph
            ASSERT([&]() {
                if(cur_idx == 0) return true;
                if(std::abs(cur_quality-quality.quality()) > EPS) {
                    LOGVAR(cur_quality);
                    LOGVAR(quality.quality());
                    return false;
                }
                return true;
            }(),"Quality of the contracted graph is not equal with quality of its corresponding uncontracted graph!");
            
            old_quality = cur_quality;
            HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
            cur_quality = louvain_pass(graph_stack[cur_idx],quality);
            HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;
            LOG("Louvain-Pass #" << iteration << " Time: " << elapsed_seconds.count() << "s");
            improvement = cur_quality - old_quality > _config.preprocessing.min_eps_improvement || max_iterations == 2;
            
            LOG("Louvain-Pass #" << iteration << " improve quality from " << old_quality << " to " << cur_quality);
            
            if(improvement) {
                cur_quality = quality.quality();
                LOG("Starting Contraction of communities...");
                start = std::chrono::high_resolution_clock::now();
                auto contraction = graph_stack[cur_idx++].contractCluster();
                end = std::chrono::high_resolution_clock::now();
                elapsed_seconds = end - start;
                LOG("Contraction Time: " << elapsed_seconds.count() << "s");
                graph_stack.push_back(std::move(contraction.first));
                mapping_stack.push_back(std::move(contraction.second));
                LOG("Current number of communities: " << graph_stack[cur_idx].numNodes());
            }
            
            LOG("");
            
        } while(improvement && iteration < max_iterations);
        
        ASSERT((mapping_stack.size() + 1) == graph_stack.size(), "Unequality between graph and mapping stack!");
        while(!mapping_stack.empty()) {
            assignClusterToNextLevelFinerGraph(graph_stack[cur_idx-1],graph_stack[cur_idx],mapping_stack[cur_idx-1]);
            graph_stack.pop_back();
            mapping_stack.pop_back();
            cur_idx--;
        }
        
        for(NodeID node : _graph.nodes()) { 
            _graph.setClusterID(node,graph_stack[0].clusterID(node));    
        }    
        
        _first_louvain_call = false;
        
        return cur_quality;
        
    }
    
    bool wasAlreadyExecuted() const {
        return !_first_louvain_call;
    }
    
    Graph getGraph() {
        return _graph;
    }
    
    ClusterID contractHypernodes(const HypernodeID hn1, const HypernodeID hn2) {
        _graph.contractHypernodes(hn1,hn2);
    }
    
    ClusterID clusterID(const HypernodeID hn) const {
        return _graph.hypernodeClusterID(hn);
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
        int iterations = 0;
        ClusterID community_limit = _config.coarsening.contraction_limit*_config.preprocessing.community_limit;
        LOGVAR(community_limit);
        
        //TODO(heuer): Think about shuffling nodes before louvain pass

        do {
            LOG("######## Starting Louvain-Pass-Iteration #" << ++iterations << " ########");
            node_moves = 0;
            for(NodeID node : g.nodes()) {
                ClusterID cur_cid = g.clusterID(node);
                EdgeWeight cur_incident_cluster_weight = 0.0L;
                ClusterID best_cid = cur_cid;
                EdgeWeight best_incident_cluster_weight = 0.0L;
                EdgeWeight best_gain = 0.0L;
              
                if(g.clusterSize(cur_cid) == 1 && g.numCommunities() <= community_limit) {
                  continue;
                }
                
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
            
            LOG("Iteration #" << iterations << ": Moving " << node_moves << " nodes to new communities.");
            
        } while(node_moves > 0 && iterations < _config.preprocessing.max_louvain_pass_iterations);
        

        return quality.quality();
    }
    
    Graph _graph;
    const Configuration& _config;
    bool _first_louvain_call;
};

}  // namespace kahypar

