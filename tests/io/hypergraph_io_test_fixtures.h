/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <cstdio>
#include <string>

#include "gmock/gmock.h"

#include "kahypar/definitions.h"
#include "kahypar/kahypar.h"
#include "kahypar/partition/coarsening/full_vertex_pair_coarsener.h"
#include "kahypar/partition/coarsening/heavy_edge_rater.h"
#include "kahypar/partition/coarsening/i_coarsener.h"
#include "kahypar/partition/partitioner.h"
#include "kahypar/partition/refinement/2way_fm_refiner.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"

using::testing::Test;

namespace kahypar {
namespace io {
class AnUnweightedHypergraphFile : public Test {
 public:
  AnUnweightedHypergraphFile() :
    _filename(""),
    _num_hyperedges(0),
    _num_hypernodes(0),
    _hypergraph_type(0),
    _control_index_vector({ 0, 2, 6, 9, 12 }),
    _control_edge_vector({ 0, 1, 0, 6, 4, 5, 4, 5, 3, 1, 2, 3 }) { }

  void SetUp() {
    _filename = "test_instances/unweighted_hypergraph.hgr";
  }

  std::string _filename;
  HyperedgeID _num_hyperedges;
  HypernodeID _num_hypernodes;
  int _hypergraph_type;
  HyperedgeIndexVector _control_index_vector;
  HyperedgeVector _control_edge_vector;
};

class AHypergraphFileWithHyperedgeWeights : public AnUnweightedHypergraphFile {
 public:
  AHypergraphFileWithHyperedgeWeights() :
    AnUnweightedHypergraphFile(),
    _control_hyperedge_weights({ 2, 3, 8, 7 }) { }

  void SetUp() {
    _filename = "test_instances/weighted_hyperedges_hypergraph.hgr";
  }

  HyperedgeWeightVector _control_hyperedge_weights;
};

class AHypergraphFileWithHypernodeWeights : public AnUnweightedHypergraphFile {
 public:
  AHypergraphFileWithHypernodeWeights() :
    AnUnweightedHypergraphFile(),
    _control_hypernode_weights({ 5, 1, 8, 7, 3, 9, 3 }) { }

  void SetUp() {
    _filename = "test_instances/weighted_hypernodes_hypergraph.hgr";
  }

  HyperedgeWeightVector _control_hypernode_weights;
};

class AHypergraphFileWithoutHyperedges : public AnUnweightedHypergraphFile {
 public:
  AHypergraphFileWithoutHyperedges() :
    AnUnweightedHypergraphFile(),
    _control_hypernode_weights({ 23, 42, 11 }) { }

  void SetUp() {
    _filename = "test_instances/hypergraph_without_hyperedges.hgr";
  }

  HypernodeWeightVector _control_hypernode_weights;
};

class AHypergraphFileWithHypernodeAndHyperedgeWeights : public AnUnweightedHypergraphFile {
 public:
  AHypergraphFileWithHypernodeAndHyperedgeWeights() :
    AnUnweightedHypergraphFile(),
    _control_hyperedge_weights({ 2, 3, 8, 7 }),
    _control_hypernode_weights({ 5, 1, 8, 7, 3, 9, 3 }) { }

  void SetUp() {
    _filename = "test_instances/weighted_hyperedges_and_hypernodes_hypergraph.hgr";
  }

  HyperedgeWeightVector _control_hyperedge_weights;
  HypernodeWeightVector _control_hypernode_weights;
};

class AnUnweightedHypergraph : public Test {
 public:
  AnUnweightedHypergraph() :
    _filename(""),
    _num_hyperedges(4),
    _num_hypernodes(7),
    _hypergraph_type(0),
    _index_vector({ 0, 2, 6, 9, 12 }),
    _edge_vector({ 0, 1, 0, 6, 4, 5, 4, 5, 3, 1, 2, 3 }),
    _written_index_vector(),
    _written_edge_vector(),
    _written_file(),
    _hypergraph(nullptr) { }

  void SetUp() {
    _filename = "test_instances/unweighted_hypergraph.hgr.out";
    _hypergraph = new Hypergraph(_num_hypernodes, _num_hyperedges, _index_vector, _edge_vector);
  }

  void TearDown() {
    std::remove(_filename.c_str());
    delete _hypergraph;
  }

  std::string _filename;
  HyperedgeID _num_hyperedges;
  HypernodeID _num_hypernodes;
  int _hypergraph_type;
  HyperedgeIndexVector _index_vector;
  HyperedgeVector _edge_vector;
  HyperedgeIndexVector _written_index_vector;
  HyperedgeVector _written_edge_vector;
  std::ifstream _written_file;
  Hypergraph* _hypergraph;
};

class AHypergraphWithHyperedgeWeights : public AnUnweightedHypergraph {
 public:
  AHypergraphWithHyperedgeWeights() :
    AnUnweightedHypergraph(),
    _hyperedge_weights({ 2, 3, 8, 7 }),
    _written_hyperedge_weights() { }

  void SetUp() {
    _filename = "test_instances/weighted_hyperedges_hypergraph.hgr.out";
    _hypergraph = new Hypergraph(_num_hypernodes, _num_hyperedges, _index_vector, _edge_vector,
                                 2, &_hyperedge_weights);
  }

  HyperedgeWeightVector _hyperedge_weights;
  HyperedgeWeightVector _written_hyperedge_weights;
};

class AHypergraphWithHypernodeWeights : public AnUnweightedHypergraph {
 public:
  AHypergraphWithHypernodeWeights() :
    AnUnweightedHypergraph(),
    _hypernode_weights({ 5, 1, 8, 7, 3, 9, 3 }),
    _written_hypernode_weights() { }

  void SetUp() {
    _filename = "test_instances/weighted_hypernodes_hypergraph.hgr.out";
    _hypergraph = new Hypergraph(_num_hypernodes, _num_hyperedges, _index_vector, _edge_vector,
                                 2, nullptr, &_hypernode_weights);
  }

  HypernodeWeightVector _hypernode_weights;
  HypernodeWeightVector _written_hypernode_weights;
};

class AHypergraphWithHypernodeAndHyperedgeWeights : public AnUnweightedHypergraph {
 public:
  AHypergraphWithHypernodeAndHyperedgeWeights() :
    AnUnweightedHypergraph(),
    _hyperedge_weights({ 2, 3, 8, 7 }),
    _hypernode_weights({ 5, 1, 8, 7, 3, 9, 3 }),
    _written_hyperedge_weights(),
    _written_hypernode_weights() { }

  void SetUp() {
    _filename = "test_instances/weighted_hyperedges_and_hypernodes_hypergraph.hgr.out";
    _hypergraph = new Hypergraph(_num_hypernodes, _num_hyperedges, _index_vector, _edge_vector,
                                 2, &_hyperedge_weights, &_hypernode_weights);
  }

  HyperedgeWeightVector _hyperedge_weights;
  HypernodeWeightVector _hypernode_weights;
  HyperedgeWeightVector _written_hyperedge_weights;
  HypernodeWeightVector _written_hypernode_weights;
};

using FirstWinsRater = HeavyEdgeRater<RatingType, FirstRatingWins>;
using FirstWinsCoarsener = FullVertexPairCoarsener<FirstWinsRater>;
using Refiner = TwoWayFMRefiner<NumberOfFruitlessMovesStopsSearch>;

class APartitionOfAHypergraph : public Test {
 public:
  APartitionOfAHypergraph() :
    _hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9,  /*sentinel*/ 12 },
                HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 }),
    _config(),
    _partitioner(),
    _coarsener(new FirstWinsCoarsener(_hypergraph, _config,  /* heaviest_node_weight */ 1)),
    _refiner(new Refiner(_hypergraph, _config)) {
    _config.partition.k = 2;
    _config.partition.rb_lower_k = 0;
    _config.partition.rb_upper_k = _config.partition.k - 1;
    _config.partition.total_graph_weight = 7;
    _config.local_search.algorithm = RefinementAlgorithm::twoway_fm;
    _config.coarsening.contraction_limit = 2;
    _config.coarsening.max_allowed_node_weight = 5;
    _config.partition.graph_filename = "APartitionOfAHypergraphTest";
    _config.partition.graph_partition_filename = "APartitionOfAHypergraphTest.hgr.part.2.KaHyPar";
    _config.partition.perfect_balance_part_weights[0] = ceil(7.0 / 2);
    _config.partition.perfect_balance_part_weights[1] = ceil(7.0 / 2);
    _config.partition.max_part_weights[0] = (1 + _config.partition.epsilon)
                                            * _config.partition.perfect_balance_part_weights[0];
    _config.partition.max_part_weights[1] = (1 + _config.partition.epsilon) *
                                            _config.partition.perfect_balance_part_weights[1];
  }

  void TearDown() {
    std::remove(_config.partition.graph_partition_filename.c_str());
  }

  Hypergraph _hypergraph;
  Configuration _config;
  Partitioner _partitioner;
  std::unique_ptr<ICoarsener> _coarsener;
  std::unique_ptr<IRefiner> _refiner;
};
}  // namespace io
}  // namespace kahypar
