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

#include "gmock/gmock.h"

#include "kahypar/definitions.h"
#include "kahypar/kahypar.h"
#include "kahypar/partition/coarsening/full_vertex_pair_coarsener.h"
#include "kahypar/partition/coarsening/heavy_edge_rater.h"
#include "kahypar/partition/coarsening/i_coarsener.h"
#include "kahypar/partition/configuration.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/partitioner.h"
#include "kahypar/partition/refinement/2way_fm_refiner.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"

using::testing::Test;
using::testing::Eq;
using::testing::DoubleEq;

namespace kahypar {
namespace metrics {
using FirstWinsRater = HeavyEdgeRater<RatingType, FirstRatingWins>;
using FirstWinsCoarsener = FullVertexPairCoarsener<FirstWinsRater>;
using Refiner = TwoWayFMRefiner<NumberOfFruitlessMovesStopsSearch>;

class AnUnPartitionedHypergraph : public Test {
 public:
  AnUnPartitionedHypergraph() :
    hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9,  /*sentinel*/ 12 },
               HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 }) { }

  Hypergraph hypergraph;
};

class TheDemoHypergraph : public AnUnPartitionedHypergraph {
 public:
  TheDemoHypergraph() :
    AnUnPartitionedHypergraph() { }
};

class APartitionedHypergraph : public Test {
 public:
  APartitionedHypergraph() :
    hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9,  /*sentinel*/ 12 },
               HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 }),
    config(),
    partitioner(),
    coarsener(new FirstWinsCoarsener(hypergraph, config,  /* heaviest_node_weight */ 1)),
    refiner(new Refiner(hypergraph, config)) {
    config.partition.k = 2;
    config.local_search.algorithm = RefinementAlgorithm::twoway_fm;
    config.coarsening.contraction_limit = 2;
    config.partition.total_graph_weight = 7;
    config.coarsening.max_allowed_node_weight = 5;
    config.partition.graph_filename = "Test";
    config.partition.graph_partition_filename = "Test.hgr.part.2.KaHyPar";
    config.partition.epsilon = 0.15;
    config.partition.perfect_balance_part_weights[0] = ceil(7.0 / 2);
    config.partition.perfect_balance_part_weights[1] = ceil(7.0 / 2);
    config.partition.max_part_weights[0] = (1 + config.partition.epsilon)
                                           * config.partition.perfect_balance_part_weights[0];
    config.partition.max_part_weights[1] = (1 + config.partition.epsilon)
                                           * config.partition.perfect_balance_part_weights[1];
    partitioner.performPartitioning(hypergraph, *coarsener, *refiner, config);
  }

  Hypergraph hypergraph;
  Configuration config;
  Partitioner partitioner;
  std::unique_ptr<ICoarsener> coarsener;
  std::unique_ptr<IRefiner> refiner;
};

class TheHyperedgeCutCalculationForInitialPartitioning : public AnUnPartitionedHypergraph {
 public:
  TheHyperedgeCutCalculationForInitialPartitioning() :
    AnUnPartitionedHypergraph(),
    config(),
    coarsener(hypergraph, config,  /* heaviest_node_weight */ 1),
    hg_to_hmetis(),
    partition() {
    config.coarsening.contraction_limit = 2;
    config.coarsening.max_allowed_node_weight = 5;
    config.partition.graph_filename = "cutCalc_test.hgr";
    config.partition.graph_partition_filename = "cutCalc_test.hgr.part.2.KaHyPar";
    config.partition.epsilon = 0.15;
    hg_to_hmetis[1] = 0;
    hg_to_hmetis[3] = 1;
    partition.push_back(1);
    partition.push_back(0);
  }

  Configuration config;
  FirstWinsCoarsener coarsener;
  std::unordered_map<HypernodeID, HypernodeID> hg_to_hmetis;
  std::vector<PartitionID> partition;
};

TEST_F(TheHyperedgeCutCalculationForInitialPartitioning, ReturnsCorrectResult) {
  coarsener.coarsen(2);
  ASSERT_THAT(hypergraph.nodeDegree(1), Eq(1));
  ASSERT_THAT(hypergraph.nodeDegree(3), Eq(1));
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(3, 1);

  ASSERT_THAT(hyperedgeCut(hypergraph, hg_to_hmetis, partition), Eq(hyperedgeCut(hypergraph)));
}

TEST_F(AnUnPartitionedHypergraph, HasHyperedgeCutZero) {
  ASSERT_THAT(hyperedgeCut(hypergraph), Eq(0));
}

TEST_F(APartitionedHypergraph, HasCorrectHyperedgeCut) {
  ASSERT_THAT(hyperedgeCut(hypergraph), Eq(2));
}

TEST_F(TheDemoHypergraph, HasAvgHyperedgeDegree3) {
  ASSERT_THAT(avgHyperedgeDegree(hypergraph), DoubleEq(3.0));
}

TEST_F(TheDemoHypergraph, HasAvgHypernodeDegree12Div7) {
  ASSERT_THAT(avgHypernodeDegree(hypergraph), DoubleEq(12.0 / 7));
}
}  // namespace metrics
}  // namespace kahypar
