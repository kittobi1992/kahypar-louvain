/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#include "gmock/gmock.h"

#include "lib/definitions.h"

using defs::Hypergraph;
using defs::HyperedgeIndexVector;
using defs::HyperedgeVector;
using defs::HyperedgeID;

namespace datastructure {

TEST(AHypergraphNeighborhood, MergesNeighborhoodsOnContraction) {
  Hypergraph hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9,  /*sentinel*/ 12 },
                        HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 });

  ASSERT_THAT(hypergraph.neighborhood().of(0),
              ::testing::ContainerEq(std::vector<int>{0,1,2,3,4}));
  ASSERT_THAT(hypergraph.neighborhood().of(3),
              ::testing::ContainerEq(std::vector<int>{0,1,3,4,6}));

  hypergraph.contract(0,3);

  ASSERT_THAT(hypergraph.neighborhood().of(0),
              ::testing::ContainerEq(std::vector<int>{0,1,2,4,6}));
}

TEST(AHypergraphNeighborhood, SupportsUncontraction) {
  Hypergraph hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9,  /*sentinel*/ 12 },
                        HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 });

  ASSERT_THAT(hypergraph.neighborhood().of(0),
              ::testing::ContainerEq(std::vector<int>{0,1,2,3,4}));
  ASSERT_THAT(hypergraph.neighborhood().of(3),
              ::testing::ContainerEq(std::vector<int>{0,1,3,4,6}));

  auto memento = hypergraph.contract(0,3);
  hypergraph.setNodePart(0,0);
  hypergraph.uncontract(memento);

  ASSERT_THAT(hypergraph.neighborhood().of(0),
              ::testing::ContainerEq(std::vector<int>{0,1,2,3,4}));
  ASSERT_THAT(hypergraph.neighborhood().of(3),
              ::testing::ContainerEq(std::vector<int>{0,1,3,4,6}));
}





TEST(AHypergraphNeighborhood, RemovesContractedHypernodeFromNeighborhoods) {
  Hypergraph hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9,  /*sentinel*/ 12 },
                        HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 });

  ASSERT_THAT(hypergraph.neighborhood().of(4),
              ::testing::ContainerEq(std::vector<int>{0,1,3,4,6}));
  ASSERT_THAT(hypergraph.neighborhood().of(3),
              ::testing::ContainerEq(std::vector<int>{0,1,3,4,6}));

  hypergraph.contract(0,3);

  ASSERT_THAT(hypergraph.neighborhood().of(4),
              ::testing::ContainerEq(std::vector<int>{0,1,4,6}));
  ASSERT_THAT(hypergraph.neighborhood().of(3),
              ::testing::ContainerEq(std::vector<int>{0,1,4,6}));
}


}
