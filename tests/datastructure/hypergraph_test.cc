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

#include <iostream>
#include <stack>
#include <tuple>

#include "gmock/gmock.h"

#include "kahypar/definitions.h"
#include "kahypar/partition/coarsening/coarsening_memento.h"
#include "kahypar/partition/coarsening/hypergraph_pruner.h"
#include "tests/datastructure/hypergraph_test_fixtures.h"

using::testing::Eq;
using::testing::ContainerEq;
using::testing::Test;

namespace kahypar {
namespace ds {
using Memento = Hypergraph::ContractionMemento;
TEST_F(AHypergraph, InitializesInternalHypergraphRepresentation) {
  ASSERT_THAT(hypergraph.currentNumNodes(), Eq(7));
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(4));
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(12));
  ASSERT_THAT(hypergraph.nodeDegree(0), Eq(2));
  ASSERT_THAT(hypergraph.nodeDegree(1), Eq(1));
  ASSERT_THAT(hypergraph.nodeDegree(2), Eq(2));
  ASSERT_THAT(hypergraph.nodeDegree(3), Eq(2));
  ASSERT_THAT(hypergraph.nodeDegree(4), Eq(2));
  ASSERT_THAT(hypergraph.nodeDegree(5), Eq(1));
  ASSERT_THAT(hypergraph.nodeDegree(6), Eq(2));

  ASSERT_THAT(hypergraph.edgeSize(0), Eq(2));
  ASSERT_THAT(hypergraph.edgeSize(1), Eq(4));
  ASSERT_THAT(hypergraph.edgeSize(2), Eq(3));
  ASSERT_THAT(hypergraph.edgeSize(3), Eq(3));
}

TEST_F(AHypergraph, ReturnsHyperNodeDegree) {
  ASSERT_THAT(hypergraph.nodeDegree(6), Eq(2));
}

TEST_F(AHypergraph, ReturnsHyperEdgeSize) {
  ASSERT_THAT(hypergraph.edgeSize(2), Eq(3));
}

TEST_F(AHypergraph, SetsAndGetsHyperNodeWeight) {
  hypergraph.setNodeWeight(0, 42);
  ASSERT_THAT(hypergraph.nodeWeight(0), Eq(42));
}

TEST_F(AHypergraph, SetsAndGetsHyperEdgeWeight) {
  hypergraph.setEdgeWeight(1, 23);
  ASSERT_THAT(hypergraph.edgeWeight(1), Eq(23));
}

TEST_F(AHypergraph, ReturnsNumberOfHypernodes) {
  ASSERT_THAT(hypergraph.currentNumNodes(), Eq(7));
}

TEST_F(AHypergraph, ReturnsNumberOfHyperedges) {
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(4));
}

TEST_F(AHypergraph, ReturnsNumberOfPins) {
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(12));
}

TEST_F(AHypergraph, DecrementsNumberOfHypernodesOnHypernodeRemoval) {
  ASSERT_THAT(hypergraph.currentNumNodes(), Eq(7));
  hypergraph.removeNode(6);
  ASSERT_THAT(hypergraph.currentNumNodes(), Eq(6));
}

TEST_F(AHypergraph, DecrementsNumberOfPinsOnHypernodeRemoval) {
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(12));
  hypergraph.removeNode(6);
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(10));
}

TEST_F(AHypergraph, DecrementsSizeOfAffectedHyperedgesOnHypernodeRemoval) {
  ASSERT_THAT(hypergraph.edgeSize(3), Eq(3));
  ASSERT_THAT(hypergraph.edgeSize(2), Eq(3));
  hypergraph.removeNode(6);
  ASSERT_THAT(hypergraph.edgeSize(3), Eq(2));
  ASSERT_THAT(hypergraph.edgeSize(2), Eq(2));
}

TEST_F(AHypergraph, InvalidatesRemovedHypernode) {
  ASSERT_THAT(hypergraph.nodeIsEnabled(6), Eq(true));
  hypergraph.removeNode(6);
  ASSERT_THAT(hypergraph.nodeIsEnabled(6), Eq(false));
}

TEST_F(AHypergraph, DecrementsNumberOfHyperedgesOnHyperedgeRemoval) {
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(4));
  hypergraph.removeEdge(2);
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(3));
}

TEST_F(AHypergraph, InvalidatesRemovedHyperedge) {
  ASSERT_THAT(hypergraph.edgeIsEnabled(2), Eq(true));
  hypergraph.removeEdge(2);
  ASSERT_THAT(hypergraph.edgeIsEnabled(2), Eq(false));
}

TEST_F(AHypergraph, DecrementsNumberOfPinsOnHyperedgeRemoval) {
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(12));
  hypergraph.removeEdge(2);
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(9));
}

TEST_F(AHypergraph, DecrementsHypernodeDegreeOfAffectedHypernodesOnHyperedgeRemoval) {
  ASSERT_THAT(hypergraph.hypernode(3).size(), Eq(2));
  ASSERT_THAT(hypergraph.hypernode(4).size(), Eq(2));
  ASSERT_THAT(hypergraph.hypernode(6).size(), Eq(2));
  hypergraph.removeEdge(2);
  ASSERT_THAT(hypergraph.hypernode(3).size(), Eq(1));
  ASSERT_THAT(hypergraph.hypernode(4).size(), Eq(1));
  ASSERT_THAT(hypergraph.hypernode(6).size(), Eq(1));
}

TEST_F(AHypergraph, InvalidatesContractedHypernode) {
  ASSERT_THAT(hypergraph.nodeIsEnabled(2), Eq(true));
  hypergraph.contract(0, 2);
  ASSERT_THAT(hypergraph.nodeIsEnabled(2), Eq(false));
}

TEST_F(AHypergraph, RelinksHyperedgesOfContractedHypernodeToRepresentative) {
  ASSERT_THAT(hypergraph.nodeDegree(0), Eq(2));
  hypergraph.contract(0, 2);
  hypergraph.contract(0, 4);
  ASSERT_THAT(hypergraph.nodeDegree(0), Eq(4));
}

TEST_F(AHypergraph, AddsHypernodeWeightOfContractedNodeToRepresentative) {
  ASSERT_THAT(hypergraph.nodeWeight(0), Eq(1));
  hypergraph.contract(0, 2);
  ASSERT_THAT(hypergraph.nodeWeight(0), Eq(2));
}

TEST_F(AHypergraph, ReducesHyperedgeSizeOfHyperedgesAffectedByContraction) {
  ASSERT_THAT(hypergraph.edgeSize(0), Eq(2));
  hypergraph.contract(0, 2);
  ASSERT_THAT(hypergraph.edgeSize(0), Eq(1));
}

TEST_F(AHypergraph, ReducesNumberOfPinsOnContraction) {
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(12));
  hypergraph.contract(3, 4);
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(10));
}

TEST_F(AHypergraph, ReducesTheNumberOfHypernodesOnContraction) {
  ASSERT_THAT(hypergraph.currentNumNodes(), Eq(7));
  hypergraph.contract(3, 4);
  ASSERT_THAT(hypergraph.currentNumNodes(), Eq(6));
}

TEST_F(AHypergraph, DoesNotRemoveParallelHyperedgesOnContraction) {
  ASSERT_THAT(hypergraph.nodeDegree(0), Eq(2));
  hypergraph.contract(5, 6);
  hypergraph.contract(0, 5);
  ASSERT_THAT(hypergraph.nodeDegree(0), Eq(4));
  ASSERT_THAT(hypergraph.edgeIsEnabled(0), Eq(true));
  ASSERT_THAT(hypergraph.edgeIsEnabled(3), Eq(true));
  ASSERT_THAT(hypergraph.edgeWeight(0), Eq(1));
  ASSERT_THAT(hypergraph.edgeWeight(3), Eq(1));
}

TEST_F(AHypergraph, DoesNotRemoveHyperedgesOfSizeOneOnContraction) {
  hypergraph.contract(0, 2);
  ASSERT_THAT(hypergraph.edgeSize(0), Eq(1));

  ASSERT_THAT(hypergraph.edgeIsEnabled(0), Eq(true));
}

TEST_F(AHypernodeIterator, StartsWithFirstHypernode) {
  ASSERT_THAT(*(hypergraph.nodes().first), Eq(0));
}

TEST_F(AHypernodeIterator, BeginsWithTheFirstValidWhenIterating) {
  hypergraph.removeNode(0);
  ASSERT_THAT(*(hypergraph.nodes().first), Eq(1));
}

TEST_F(AHypernodeIterator, SkipsInvalidHypernodesWhenForwardIterating) {
  hypergraph.removeNode(1);
  hypergraph.removeNode(2);
  auto begin = hypergraph.nodes().first;
  ++begin;
  ASSERT_THAT(*begin, Eq(3));
}

TEST_F(AHyperedgeIterator, StartsWithFirstHyperedge) {
  ASSERT_THAT(*(hypergraph.edges().first), Eq(0));
}

TEST_F(AHyperedgeIterator, StartsWithTheFirstValidHyperedge) {
  hypergraph.removeEdge(0);
  ASSERT_THAT(*(hypergraph.edges().first), Eq(1));
}

TEST_F(AHyperedgeIterator, SkipsInvalidHyperedgesWhenForwardIterating) {
  hypergraph.removeEdge(1);
  hypergraph.removeEdge(2);
  auto begin = hypergraph.edges().first;
  ++begin;
  ASSERT_THAT(*begin, Eq(3));
}

TEST_F(AHypergraphMacro, IteratesOverAllHypernodes) {
  HypernodeID hypernode_count = 0;

  for (const HypernodeID hn : hypergraph.nodes()) {
    ASSERT_THAT(hn, Eq(hypernode_count));
    ++hypernode_count;
  }
  ASSERT_THAT(hypernode_count, Eq(7));
}

TEST_F(AHypergraphMacro, IteratesOverAllHyperedges) {
  HyperedgeID hyperedge_count = 0;

  for (const HyperedgeID he : hypergraph.edges()) {
    ASSERT_THAT(he, Eq(hyperedge_count));
    ++hyperedge_count;
  }
  ASSERT_THAT(hyperedge_count, Eq(4));
}

TEST_F(AHypergraphMacro, IteratesOverAllIncidentHyperedges) {
  int i = 0;

  for (const HyperedgeID he : hypergraph.incidentEdges(6)) {
    ASSERT_THAT(he, Eq(*(hypergraph._incidence_array.begin() +
                         hypergraph.hypernode(6).firstEntry() + i)));
    ++i;
  }
}

TEST_F(AHypergraphMacro, IteratesOverAllPinsOfAHyperedge) {
  int i = 0;
  for (const HypernodeID pin : hypergraph.pins(2)) {
    ASSERT_THAT(pin, Eq(*(hypergraph._incidence_array.begin() +
                          hypergraph.hyperedge(2).firstEntry() + i)));
    ++i;
  }
}

TEST_F(AContractionMemento, StoresOldStateOfInvolvedHypernodes) {
  HypernodeID u_id = 4;
  HypernodeID u_offset = hypergraph.hypernode(u_id).firstEntry();
  HypernodeID u_size = hypergraph.hypernode(u_id).size();
  HypernodeID v_id = 6;

  Memento memento = hypergraph.contract(4, 6);

  ASSERT_THAT(memento.u, Eq(u_id));
  ASSERT_THAT(memento.u_first_entry, Eq(u_offset));
  ASSERT_THAT(memento.u_size, Eq(u_size));
  ASSERT_THAT(memento.v, Eq(v_id));
}

TEST_F(AnUncontractionOperation, ReEnablesTheInvalidatedHypernode) {
  Memento memento = hypergraph.contract(4, 6);

  ASSERT_THAT(hypergraph.nodeIsEnabled(6), Eq(false));

  hypergraph.uncontract(memento);

  ASSERT_THAT(hypergraph.nodeIsEnabled(6), Eq(true));
}

TEST_F(AnUncontractionOperation, ResetsWeightOfRepresentative) {
  ASSERT_THAT(hypergraph.nodeWeight(4), Eq(1));
  Memento memento = hypergraph.contract(4, 6);
  ASSERT_THAT(hypergraph.nodeWeight(4), Eq(2));

  hypergraph.uncontract(memento);

  ASSERT_THAT(hypergraph.nodeWeight(4), Eq(1));
}

TEST_F(AnUncontractionOperation, DisconnectsHyperedgesAddedToRepresenativeDuringContraction) {
  ASSERT_THAT(hypergraph.nodeDegree(4), Eq(2));
  Memento memento = hypergraph.contract(4, 6);
  ASSERT_THAT(hypergraph.nodeDegree(4), Eq(3));

  hypergraph.uncontract(memento);

  ASSERT_THAT(hypergraph.nodeDegree(4), Eq(2));
}

TEST_F(AnUncontractionOperation, DeletesIncidenceInfoAddedDuringContraction) {
  ASSERT_THAT(hypergraph._incidence_array.size(), Eq(24));
  Memento memento = hypergraph.contract(4, 6);
  ASSERT_THAT(hypergraph._incidence_array.size(), Eq(27));

  hypergraph.uncontract(memento);
  ASSERT_THAT(hypergraph._incidence_array.size(), Eq(24));
}

TEST_F(AnUncontractionOperation, RestoresIncidenceInfoForHyperedgesAddedToRepresentative) {
  ASSERT_THAT(std::count(hypergraph.pins(3).first, hypergraph.pins(3).second, 6), Eq(1));
  ASSERT_THAT(std::count(hypergraph.pins(2).first, hypergraph.pins(2).second, 6), Eq(1));
  Memento memento = hypergraph.contract(4, 6);
  ASSERT_THAT(std::count(hypergraph.pins(3).first, hypergraph.pins(3).second, 6), Eq(0));
  ASSERT_THAT(std::count(hypergraph.pins(2).first, hypergraph.pins(2).second, 6), Eq(0));

  hypergraph.uncontract(memento);

  ASSERT_THAT(std::count(hypergraph.pins(3).first, hypergraph.pins(3).second, 6), Eq(1));
  ASSERT_THAT(std::count(hypergraph.pins(2).first, hypergraph.pins(2).second, 6), Eq(1));
}

TEST_F(AnUncontractionOperation, RestoresIncidenceInfoForHyperedgesAlredyExistingAtRepresentative) {
  ASSERT_THAT(std::count(hypergraph.pins(2).first, hypergraph.pins(2).second, 4), Eq(1));
  ASSERT_THAT(std::count(hypergraph.pins(1).first, hypergraph.pins(1).second, 4), Eq(1));
  Memento memento = hypergraph.contract(3, 4);
  ASSERT_THAT(std::count(hypergraph.pins(2).first, hypergraph.pins(2).second, 4), Eq(0));
  ASSERT_THAT(std::count(hypergraph.pins(1).first, hypergraph.pins(1).second, 4), Eq(0));

  hypergraph.uncontract(memento);

  ASSERT_THAT(std::count(hypergraph.pins(2).first, hypergraph.pins(2).second, 4), Eq(1));
  ASSERT_THAT(std::count(hypergraph.pins(1).first, hypergraph.pins(1).second, 4), Eq(1));
}

TEST_F(AnUncontractionOperation, RestoresNumberOfPinsOnUncontraction) {
  hypergraph.uncontract(hypergraph.contract(3, 4));
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(12));
}

TEST_F(AnUncontractionOperation, RestoresHyperedgeSizeOfHyperedgesAffectedByContraction) {
  hypergraph.uncontract(hypergraph.contract(0, 2));
  ASSERT_THAT(hypergraph.edgeSize(0), Eq(2));
}

TEST_F(AnUncontractedHypergraph, EqualsTheInitialHypergraphBeforeContraction) {
  std::vector<std::pair<HypernodeID, HypernodeID> > contractions { { 4, 6 }, { 3, 4 }, { 0, 2 },
                                                                   { 0, 1 }, { 0, 5 }, { 0, 3 } };
  std::stack<CoarseningMemento> contraction_history;
  HypergraphPruner hypergraph_pruner(modified_hypergraph.initialNumNodes());
  for (const auto& contraction : contractions) {
    contraction_history.emplace(modified_hypergraph.contract(contraction.first,
                                                             contraction.second));
    hypergraph_pruner.removeSingleNodeHyperedges(modified_hypergraph, contraction_history.top());
  }

  ASSERT_THAT(modified_hypergraph.nodeWeight(0), Eq(7));
  modified_hypergraph.setNodePart(0, 0);

  while (!contraction_history.empty()) {
    hypergraph_pruner.restoreSingleNodeHyperedges(modified_hypergraph,
                                                  contraction_history.top());
    modified_hypergraph.uncontract(contraction_history.top().contraction_memento);
    contraction_history.pop();
  }

  ASSERT_THAT(verifyEquivalenceWithoutPartitionInfo(hypergraph, modified_hypergraph), Eq(true));
}

TEST_F(AHypergraph, ReturnsInitialNumberOfHypernodesAfterHypergraphModification) {
  ASSERT_THAT(hypergraph.initialNumNodes(), Eq(7));
  hypergraph.removeNode(6);
  ASSERT_THAT(hypergraph.initialNumNodes(), Eq(7));
}

TEST_F(AHypergraph, ReturnsInitialNumberOfPinsAfterHypergraphModification) {
  ASSERT_THAT(hypergraph.initialNumPins(), Eq(12));
  hypergraph.removeNode(6);
  ASSERT_THAT(hypergraph.initialNumPins(), Eq(12));
}

TEST_F(AHypergraph, ReturnsInitialNumberHyperedgesAfterHypergraphModification) {
  ASSERT_THAT(hypergraph.initialNumEdges(), Eq(4));
  hypergraph.removeEdge(2);
  ASSERT_THAT(hypergraph.initialNumEdges(), Eq(4));
}

TEST_F(AnUncontractionOperation, UpdatesPartitionIndexOfUncontractedNode) {
  ASSERT_THAT(hypergraph.partID(2), Eq(0));

  Memento memento = hypergraph.contract(0, 2);
  hypergraph.changeNodePart(0, 0, 1);
  hypergraph.uncontract(memento);

  ASSERT_THAT(hypergraph.partID(2), Eq(1));
}

TEST(AnUnconnectedHypernode, IsNotRemovedTogetherWithLastEdgeIfFlagIsFalse) {
  Hypergraph hypergraph(1, 1, HyperedgeIndexVector { 0,  /*sentinel*/ 1 },
                        HyperedgeVector { 0 });

  hypergraph.removeEdge(0);
  ASSERT_THAT(hypergraph.nodeIsEnabled(0), Eq(true));
}

TEST_F(AHypergraph, ReducesPinCountOfAffectedHEsOnContraction) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 0);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 0);
  hypergraph.setNodePart(6, 0);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 0), Eq(4));
  ASSERT_THAT(hypergraph.pinCountInPart(2, 0), Eq(3));
  hypergraph.contract(3, 4);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 0), Eq(3));
  ASSERT_THAT(hypergraph.pinCountInPart(2, 0), Eq(2));
}

TEST_F(AHypergraph, IncreasesPinCountOfAffectedHEsOnUnContraction) {
  hypergraph.setNodePart(0, 1);
  hypergraph.setNodePart(1, 1);
  hypergraph.setNodePart(2, 1);
  hypergraph.setNodePart(3, 1);
  hypergraph.setNodePart(4, 1);
  hypergraph.setNodePart(5, 1);
  hypergraph.setNodePart(6, 1);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(4));
  ASSERT_THAT(hypergraph.pinCountInPart(2, 1), Eq(3));
  Memento memento = hypergraph.contract(3, 4);
  hypergraph.changeNodePart(3, 1, 0);
  hypergraph.uncontract(memento);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(2));
  ASSERT_THAT(hypergraph.pinCountInPart(1, 0), Eq(2));
  ASSERT_THAT(hypergraph.pinCountInPart(2, 1), Eq(1));
  ASSERT_THAT(hypergraph.pinCountInPart(2, 0), Eq(2));
}

TEST_F(APartitionedHypergraph, StoresPartitionPinCountsForHyperedges) {
  ASSERT_THAT(hypergraph.pinCountInPart(0, 0), Eq(1));
  ASSERT_THAT(hypergraph.pinCountInPart(0, 1), Eq(1));
  ASSERT_THAT(hypergraph.pinCountInPart(1, 0), Eq(4));
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(0));
  ASSERT_THAT(hypergraph.pinCountInPart(2, 0), Eq(2));
  ASSERT_THAT(hypergraph.pinCountInPart(2, 1), Eq(1));
  ASSERT_THAT(hypergraph.pinCountInPart(3, 0), Eq(0));
  ASSERT_THAT(hypergraph.pinCountInPart(3, 1), Eq(3));
}

TEST_F(AHypergraph, InvalidatesPartitionPinCountsOnHyperedgeRemoval) {
  hypergraph.setNodePart(0, 1);
  hypergraph.setNodePart(1, 1);
  hypergraph.setNodePart(3, 1);
  hypergraph.setNodePart(4, 1);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(4));

  hypergraph.removeEdge(1);

  for (PartitionID part = 0; part < hypergraph._k; ++part) {
    // bypass pinCountInPart because of assertions
    const HypernodeID num_pins = hypergraph._pins_in_part[1 * hypergraph._k + part];
    ASSERT_THAT(num_pins, Eq(hypergraph.kInvalidCount));
  }
}

TEST_F(AHypergraph, RestoresInvalidatedPartitionPinCountsOnHyperedgeRestore) {
  hypergraph.setNodePart(0, 1);
  hypergraph.setNodePart(1, 1);
  hypergraph.setNodePart(3, 1);
  hypergraph.setNodePart(4, 1);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(4));
  hypergraph.removeEdge(1);
  hypergraph.restoreEdge(1);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(4));
}


TEST_F(AHypergraph, MaintainsPartitionPinCountsForRestoredSingleNodeOrParallelHEs) {
  hypergraph.removeEdge(1);

  // pseudo initial partitioning
  hypergraph.setNodePart(0, 1);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(3, 1);
  hypergraph.setNodePart(4, 1);

  hypergraph.restoreEdge(1);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(3));
  ASSERT_THAT(hypergraph.pinCountInPart(1, 0), Eq(1));
}

TEST_F(AHypergraph, DecrementsPartitionPinCountOnHypernodeRemoval) {
  hypergraph.setNodePart(0, 1);
  hypergraph.setNodePart(1, 1);
  hypergraph.setNodePart(3, 1);
  hypergraph.setNodePart(4, 1);
  hypergraph.setNodePart(6, 1);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(4));
  ASSERT_THAT(hypergraph.pinCountInPart(2, 1), Eq(3));
  hypergraph.removeNode(3);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(3));
  ASSERT_THAT(hypergraph.pinCountInPart(2, 1), Eq(2));
}

TEST_F(AHypergraph, UpdatesPartitionPinCountsIfANodeChangesPartition) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 0), Eq(4));
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(0));
  hypergraph.changeNodePart(1, 0, 1);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 0), Eq(3));
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(1));
}

TEST_F(AHypergraph, CalculatesPinCountsOfAHyperedge) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 0);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 0);
  hypergraph.setNodePart(6, 0);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 0), Eq(4));
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(0));
  hypergraph.changeNodePart(1, 0, 1);
  hypergraph.changeNodePart(4, 0, 1);
  ASSERT_THAT(hypergraph.pinCountInPart(1, 0), Eq(2));
  ASSERT_THAT(hypergraph.pinCountInPart(1, 1), Eq(2));
}

TEST_F(AnUnPartitionedHypergraph, HasAllNodesInInvalidPartition) {
  for (const HypernodeID hn : hypergraph.nodes()) {
    ASSERT_THAT(hypergraph.partID(hn), Eq(Hypergraph::kInvalidPartition));
  }
}

TEST_F(AHypergraph, MaintainsAConnectivitySetForEachHyperedge) {
  ASSERT_THAT(hypergraph.connectivity(0), Eq(0));
  ASSERT_THAT(hypergraph.connectivity(1), Eq(0));
  ASSERT_THAT(hypergraph.connectivity(2), Eq(0));
  ASSERT_THAT(hypergraph.connectivity(3), Eq(0));
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 1);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 1);
  hypergraph.setNodePart(6, 1);
  ASSERT_THAT(hypergraph.connectivity(0), Eq(2));
  ASSERT_THAT(hypergraph.connectivity(1), Eq(1));
  ASSERT_THAT(hypergraph.connectivity(2), Eq(2));
  ASSERT_THAT(hypergraph.connectivity(3), Eq(1));
}

TEST_F(AHypergraph, RemovesPartFromConnectivitySetIfHEDoesNotConnectThatPart) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 1);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 1);
  hypergraph.setNodePart(6, 1);
  hypergraph.initializeNumCutHyperedges();
  ASSERT_THAT(hypergraph.connectivity(0), Eq(2));

  hypergraph.changeNodePart(2, 1, 0);

  ASSERT_THAT(hypergraph.connectivity(0), Eq(1));
}

TEST_F(AHypergraph, AddsPartToConnectivitySetIfHEConnectsThatPart) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 0);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 1);
  hypergraph.setNodePart(6, 1);
  hypergraph.initializeNumCutHyperedges();
  ASSERT_THAT(hypergraph.connectivity(0), Eq(1));

  hypergraph.changeNodePart(2, 0, 1);

  ASSERT_THAT(hypergraph.connectivity(0), Eq(2));
}


TEST_F(AHypergraph, AllowsIterationOverConnectivitySetOfAHyperege) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(2, 1);
  ASSERT_THAT(hypergraph.connectivity(0), Eq(2));

  int part_id = 0;
  for (const PartitionID part : hypergraph.connectivitySet(0)) {
    ASSERT_THAT(part, Eq(part_id));
    ++part_id;
  }
}

TEST(ConnectivitySets, AreCleardWhenSingleNodeHyperedgesAreRemoved) {
  Hypergraph hypergraph(1, 1, HyperedgeIndexVector { 0,  /*sentinel*/ 1 },
                        HyperedgeVector { 0 });
  hypergraph.setNodePart(0, 0);
  ASSERT_THAT(hypergraph.connectivity(0), Eq(1));
  ASSERT_THAT(*hypergraph.connectivitySet(0).begin(), Eq(0));

  hypergraph.removeEdge(0);
  hypergraph.changeNodePart(0, 0, 1);
  hypergraph.restoreEdge(0);

  ASSERT_THAT(hypergraph.connectivity(0), Eq(1));
  ASSERT_THAT(*hypergraph.connectivitySet(0).begin(), Eq(1));
}

TEST_F(AHypergraph, MaintainsCorrectPartSizesDuringUncontraction) {
  std::stack<Memento> mementos;
  mementos.push(hypergraph.contract(0, 1));
  mementos.push(hypergraph.contract(0, 3));
  mementos.push(hypergraph.contract(0, 4));
  mementos.push(hypergraph.contract(2, 5));
  mementos.push(hypergraph.contract(2, 6));
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(2, 1);
  hypergraph.initializeNumCutHyperedges();
  ASSERT_THAT(hypergraph.partSize(0), Eq(1));
  ASSERT_THAT(hypergraph.partSize(1), Eq(1));

  hypergraph.uncontract(mementos.top());
  mementos.pop();
  ASSERT_THAT(hypergraph.partSize(1), Eq(2));
  hypergraph.uncontract(mementos.top());
  mementos.pop();
  ASSERT_THAT(hypergraph.partSize(1), Eq(3));
  hypergraph.uncontract(mementos.top());
  mementos.pop();
  ASSERT_THAT(hypergraph.partSize(0), Eq(2));
  hypergraph.uncontract(mementos.top());
  mementos.pop();
  ASSERT_THAT(hypergraph.partSize(0), Eq(3));
  hypergraph.uncontract(mementos.top());
  mementos.pop();
  ASSERT_THAT(hypergraph.partSize(0), Eq(4));
}

TEST_F(AHypergraph, MaintainsItsTotalWeight) {
  ASSERT_THAT(hypergraph.totalWeight(), Eq(7));
}

TEST_F(APartitionedHypergraph, CanBeDecomposedIntoHypergraphs) {
  auto extr_part0 = extractPartAsUnpartitionedHypergraphForBisection(hypergraph, 0);
  auto extr_part1 = extractPartAsUnpartitionedHypergraphForBisection(hypergraph, 1);
  Hypergraph& part0_hypergraph = *extr_part0.first;
  Hypergraph& part1_hypergraph = *extr_part1.first;

  ASSERT_THAT(part0_hypergraph.initialNumNodes(), Eq(4));
  ASSERT_THAT(part1_hypergraph.initialNumNodes(), Eq(3));

  ASSERT_THAT(part0_hypergraph.currentNumEdges(), Eq(1));
  ASSERT_THAT(part1_hypergraph.currentNumEdges(), Eq(1));

  ASSERT_THAT(part0_hypergraph.edgeSize(0), Eq(4));
  ASSERT_THAT(part1_hypergraph.edgeSize(0), Eq(3));

  const std::vector<HypernodeID>& mapping_0 = extr_part0.second;
  const std::vector<HypernodeID>& mapping_1 = extr_part1.second;

  ASSERT_THAT(mapping_0, ContainerEq(std::vector<HypernodeID>{ 0, 1, 3, 4 }));
  ASSERT_THAT(mapping_1, ContainerEq(std::vector<HypernodeID>{ 2, 5, 6 }));
}

TEST_F(AHypergraph, WithOnePartitionEqualsTheExtractedHypergraphExceptForPartitionRelatedInfos) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 0);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 0);
  hypergraph.setNodePart(6, 0);
  auto extr_part0 = extractPartAsUnpartitionedHypergraphForBisection(hypergraph, 0);
  ASSERT_THAT(verifyEquivalenceWithoutPartitionInfo(hypergraph, *extr_part0.first), Eq(true));
}

TEST_F(AHypergraph, ExtractedFromAPartitionedHypergraphHasInitializedPartitionInformation) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 0);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 0);
  hypergraph.setNodePart(6, 0);
  auto extr_part0 = extractPartAsUnpartitionedHypergraphForBisection(hypergraph, 0);

  ASSERT_THAT(extr_part0.first->_part_info.size(), Eq(2));
  ASSERT_THAT(extr_part0.first->_pins_in_part.size(), Eq(8));

  for (const HyperedgeID he: extr_part0.first->edges()) {
    ASSERT_THAT(extr_part0.first->connectivity(he), Eq(0));
  }

  for (const HypernodeID hn : extr_part0.first->nodes()) {
    ASSERT_THAT(extr_part0.first->partID(hn), Eq(-1));
  }
}


TEST_F(AHypergraph, ExtractedFromAPartitionedHypergraphHasCorrectNumberOfHyperedges) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 1);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 1);
  hypergraph.setNodePart(6, 0);
  auto extr_part0 = extractPartAsUnpartitionedHypergraphForBisection(hypergraph, 0);
  ASSERT_THAT(extr_part0.first->currentNumEdges(), Eq(2));
}

TEST_F(AHypergraph, ExtractedFromAPartitionedHypergraphHasCorrectNumberOfHypernodes) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 1);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 1);
  hypergraph.setNodePart(6, 0);
  auto extr_part0 = extractPartAsUnpartitionedHypergraphForBisection(hypergraph, 0);
  ASSERT_THAT(extr_part0.first->initialNumNodes(), Eq(5));
}

TEST_F(AHypergraph, CanBeDecomposedIntoHypergraphsUsingCutNetSplitting) {
  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(1, 0);
  hypergraph.setNodePart(2, 1);
  hypergraph.setNodePart(3, 0);
  hypergraph.setNodePart(4, 0);
  hypergraph.setNodePart(5, 1);
  hypergraph.setNodePart(6, 1);
  auto extr_part0 = extractPartAsUnpartitionedHypergraphForBisection(hypergraph, 0, true);
  auto extr_part1 = extractPartAsUnpartitionedHypergraphForBisection(hypergraph, 1, true);
  Hypergraph& part0_hypergraph = *extr_part0.first;
  Hypergraph& part1_hypergraph = *extr_part1.first;

  ASSERT_THAT(part0_hypergraph.initialNumNodes(), Eq(4));
  ASSERT_THAT(part1_hypergraph.initialNumNodes(), Eq(3));

  ASSERT_THAT(part0_hypergraph.currentNumEdges(), Eq(2));
  ASSERT_THAT(part1_hypergraph.currentNumEdges(), Eq(1));

  ASSERT_THAT(part0_hypergraph.edgeSize(0), Eq(4));
  ASSERT_THAT(part0_hypergraph.edgeSize(1), Eq(2));
  ASSERT_THAT(part1_hypergraph.edgeSize(0), Eq(3));

  const std::vector<HypernodeID>& mapping_0 = extr_part0.second;
  const std::vector<HypernodeID>& mapping_1 = extr_part1.second;

  ASSERT_THAT(mapping_0, ContainerEq(std::vector<HypernodeID>{ 0, 1, 3, 4 }));
  ASSERT_THAT(mapping_1, ContainerEq(std::vector<HypernodeID>{ 2, 5, 6 }));
}

TEST_F(AHypergraph, CreatedViaReindexingIsACopyOfTheOriginalHypergraph) {
  verifyEquivalenceWithoutPartitionInfo(hypergraph, *reindex(hypergraph).first);
}

TEST_F(AHypergraph, WithContractedHypernodesCanBeReindexed) {
  hypergraph.contract(1, 4);
  hypergraph.contract(0, 2);
  hypergraph.removeEdge(1);

  auto reindexed = reindex(hypergraph);

  ASSERT_THAT(reindexed.first->initialNumNodes(), Eq(5));
  ASSERT_THAT(reindexed.first->currentNumEdges(), Eq(3));
  ASSERT_THAT(reindexed.second.size(), Eq(5));
  ASSERT_THAT(reindexed.second, ContainerEq(std::vector<HypernodeID>{ 0, 1, 3, 5, 6 }));
}

TEST_F(APartitionedHypergraph, CanBeResetToUnpartitionedState) {
  hypergraph.resetPartitioning();
  ASSERT_THAT(verifyEquivalenceWithPartitionInfo(hypergraph, original_hypergraph), Eq(true));
}


TEST_F(APartitionedHypergraph, IdentifiesBorderHypernodes) {
  hypergraph.initializeNumCutHyperedges();
  ASSERT_THAT(hypergraph.isBorderNode(0), Eq(true));
  ASSERT_THAT(hypergraph.isBorderNode(1), Eq(false));
  ASSERT_THAT(hypergraph.isBorderNode(2), Eq(true));
  ASSERT_THAT(hypergraph.isBorderNode(3), Eq(true));
  ASSERT_THAT(hypergraph.isBorderNode(4), Eq(true));
  ASSERT_THAT(hypergraph.isBorderNode(5), Eq(false));
  ASSERT_THAT(hypergraph.isBorderNode(6), Eq(true));
}

TEST_F(AHypergraph, SupportsIsolationOfHypernodes) {
  ASSERT_THAT(hypergraph.nodeDegree(0), Eq(2));
  ASSERT_THAT(hypergraph.edgeSize(0), Eq(2));
  ASSERT_THAT(hypergraph.edgeSize(1), Eq(4));
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(12));

  const HyperedgeID old_degree = hypergraph.isolateNode(0);

  ASSERT_THAT(old_degree, Eq(2));
  ASSERT_THAT(hypergraph.nodeDegree(0), Eq(0));
  ASSERT_THAT(hypergraph.edgeSize(0), Eq(1));
  ASSERT_THAT(hypergraph.edgeSize(1), Eq(3));
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(10));
  ASSERT_THAT(std::find(hypergraph.pins(0).first, hypergraph.pins(0).second, 0) ==
              hypergraph.pins(0).second, Eq(true));
  ASSERT_THAT(std::find(hypergraph.pins(1).first, hypergraph.pins(1).second, 0) ==
              hypergraph.pins(1).second, Eq(true));
}

TEST_F(AHypergraph, RemovesEmptyHyperedgesOnHypernodeIsolation) {
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(4));
  ASSERT_THAT(!hypergraph.hyperedge(0).isDisabled(), Eq(true));

  hypergraph.isolateNode(0);
  hypergraph.isolateNode(2);

  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(3));
  ASSERT_THAT(hypergraph.hyperedge(0).isDisabled(), Eq(true));
}

TEST_F(AHypergraph, RestoresRemovedEmptyHyperedgesOnRestoreOfIsolatedHypernodes) {
  ASSERT_THAT(*hypergraph.pins(0).first, Eq(0));
  ASSERT_THAT(*(hypergraph.pins(0).first + 1), Eq(2));
  const HyperedgeID old_degree_0 = hypergraph.isolateNode(0);
  const HyperedgeID old_degree_2 = hypergraph.isolateNode(2);
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(3));
  ASSERT_THAT(hypergraph.hyperedge(0).isDisabled(), Eq(true));

  hypergraph.setNodePart(0, 0);
  hypergraph.setNodePart(2, 0);
  // reverse order!
  hypergraph.reconnectIsolatedNode(2, old_degree_2);
  hypergraph.reconnectIsolatedNode(0, old_degree_0);

  ASSERT_THAT(hypergraph.edgeSize(0), Eq(2));
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(4));
  ASSERT_THAT(hypergraph.hyperedge(0).isDisabled(), Eq(false));
  ASSERT_THAT(*hypergraph.pins(0).first, Eq(2));
  ASSERT_THAT(*(hypergraph.pins(0).first + 1), Eq(0));
}

TEST_F(AHypergraph, SupportsRestoreOfIsolatedHypernodes) {
  const HyperedgeID old_degree = hypergraph.isolateNode(0);
  ASSERT_THAT(hypergraph.nodeDegree(0), Eq(0));
  ASSERT_THAT(hypergraph.edgeSize(0), Eq(1));
  ASSERT_THAT(hypergraph.edgeSize(1), Eq(3));
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(10));

  hypergraph.setNodePart(0, 0);
  hypergraph.reconnectIsolatedNode(0, old_degree);

  ASSERT_THAT(hypergraph.nodeDegree(0), Eq(2));
  ASSERT_THAT(hypergraph.edgeSize(0), Eq(2));
  ASSERT_THAT(hypergraph.edgeSize(1), Eq(4));
  ASSERT_THAT(hypergraph.currentNumPins(), Eq(12));
  ASSERT_THAT(std::find(hypergraph.pins(0).first, hypergraph.pins(0).second, 0) !=
              hypergraph.pins(0).second, Eq(true));
  ASSERT_THAT(std::find(hypergraph.pins(1).first, hypergraph.pins(1).second, 0) !=
              hypergraph.pins(1).second, Eq(true));
}
}  // namespace ds
}  // namespace kahypar
