/***************************************************************************
 *  Copyright (C) 2015 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#include "gmock/gmock.h"

#include "lib/definitions.h"
#include "partition/refinement/KWayFMRefiner.h"
#include "partition/refinement/policies/FMStopPolicies.h"

using::testing::Test;
using::testing::Eq;

using defs::Hypergraph;
using defs::HyperedgeIndexVector;
using defs::HyperedgeVector;
using defs::HyperedgeWeight;
using defs::HypernodeID;

namespace partition {
using KWayFMRefinerSimpleStopping = KWayFMRefiner<NumberOfFruitlessMovesStopsSearch>;

class AKwayFMRefiner : public Test {
 public:
  AKwayFMRefiner() :
    config(),
    hypergraph(new Hypergraph(2, 2, HyperedgeIndexVector { 0, 2,  /*sentinel*/ 3 },
                              HyperedgeVector { 0, 1, 0 }, 2)),
    refiner() {
    config.fm_local_search.max_number_of_fruitless_moves = 50;
    config.partition.total_graph_weight = 2;
    config.partition.k = 2;
    config.partition.rb_lower_k = 0;
    config.partition.rb_upper_k = config.partition.k - 1;
    config.partition.epsilon = 1.0;
    config.partition.perfect_balance_part_weights[0] = ceil(config.partition.total_graph_weight /
                                                            static_cast<double>(config.partition.k));
    config.partition.perfect_balance_part_weights[1] = ceil(config.partition.total_graph_weight /
                                                            static_cast<double>(config.partition.k));

    config.partition.max_part_weights[0] = (1 + config.partition.epsilon)
                                           * config.partition.perfect_balance_part_weights[0];
    config.partition.max_part_weights[1] = config.partition.max_part_weights[0];

    hypergraph->setNodePart(0, 0);
    hypergraph->setNodePart(1, 1);
    hypergraph->initializeNumCutHyperedges();

    refiner = std::make_unique<KWayFMRefinerSimpleStopping>(*hypergraph, config);
#ifdef USE_BUCKET_PQ
    refiner->initialize(100);
#else
    refiner->initialize();
#endif
  }

  Configuration config;
  std::unique_ptr<Hypergraph> hypergraph;
  std::unique_ptr<KWayFMRefinerSimpleStopping> refiner;
};

using AKwayFMRefinerDeathTest = AKwayFMRefiner;

TEST_F(AKwayFMRefinerDeathTest, ConsidersSingleNodeHEsDuringInitialGainComputation) {
  ASSERT_DEBUG_DEATH(refiner->insertHNintoPQ(0, 10), ".*");
}

TEST_F(AKwayFMRefinerDeathTest, ConsidersSingleNodeHEsDuringInducedGainComputation) {
  ASSERT_DEBUG_DEATH(refiner->gainInducedByHypergraph(0, 0), ".*");
}

TEST_F(AKwayFMRefiner, KnowsIfAHyperedgeIsFullyActive) {
  hypergraph.reset(new Hypergraph(3, 1, HyperedgeIndexVector { 0,  /*sentinel*/ 3 },
                                  HyperedgeVector { 0, 1, 2 }, 2));
  hypergraph->setNodePart(0, 0);
  hypergraph->setNodePart(1, 0);
  hypergraph->setNodePart(2, 0);
  hypergraph->initializeNumCutHyperedges();

  refiner = std::make_unique<KWayFMRefinerSimpleStopping>(*hypergraph, config);
  refiner->initialize(100);

  refiner->_hg.activate(0);
  hypergraph->changeNodePart(0, 0, 1);
  refiner->_hg.mark(0);

  refiner->fullUpdate(0, 0, 1, 0, 42);
  ASSERT_THAT(refiner->_he_fully_active[0], Eq(true));
}
}  // namespace partition
