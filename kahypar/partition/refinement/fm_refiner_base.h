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

#include <limits>
#include <vector>

#include "kahypar/datastructure/kway_priority_queue.h"
#include "kahypar/definitions.h"
#include "kahypar/partition/configuration.h"

namespace kahypar {
static const bool dbg_refinement_fm_border_node_check = false;
static const bool dbg_refinement_kway_fm_move = false;

struct RollbackInfo {
  HypernodeID hn;
  PartitionID from_part;
  PartitionID to_part;
};

template <typename RollbackElement = Mandatory>
class FMRefinerBase {
 protected:
  static constexpr HypernodeID kInvalidHN = std::numeric_limits<HypernodeID>::max();
  static constexpr Gain kInvalidGain = std::numeric_limits<Gain>::min();
  static constexpr HyperedgeWeight kInvalidDecrease = std::numeric_limits<PartitionID>::min();

  enum HEState {
    free = std::numeric_limits<PartitionID>::max() - 1,
    locked = std::numeric_limits<PartitionID>::max(),
  };

  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain,
                                                 std::numeric_limits<Gain> >;


  FMRefinerBase(Hypergraph& hypergraph, const Configuration& config) :
    _hg(hypergraph),
    _config(config),
    _pq(config.partition.k),
    _performed_moves(),
    _hns_to_activate() {
    _performed_moves.reserve(_hg.initialNumNodes());
    _hns_to_activate.reserve(_hg.initialNumNodes());
  }

  ~FMRefinerBase() { }

  FMRefinerBase(const FMRefinerBase&) = delete;
  FMRefinerBase& operator= (const FMRefinerBase&) = delete;

  FMRefinerBase(FMRefinerBase&&) = delete;
  FMRefinerBase& operator= (FMRefinerBase&&) = delete;

  bool hypernodeIsConnectedToPart(const HypernodeID pin, const PartitionID part) const {
    for (const HyperedgeID he : _hg.incidentEdges(pin)) {
      if (_hg.pinCountInPart(he, part) > 0) {
        return true;
      }
    }
    return false;
  }

  bool moveIsFeasible(const HypernodeID max_gain_node, const PartitionID from_part,
                      const PartitionID to_part) const {
    ASSERT(_config.partition.mode == Mode::direct_kway,
           "Method should only be called in direct partitioning");
    return (_hg.partWeight(to_part) + _hg.nodeWeight(max_gain_node)
            <= _config.partition.max_part_weights[0]) && (_hg.partSize(from_part) - 1 != 0);
  }

  void moveHypernode(const HypernodeID hn, const PartitionID from_part,
                     const PartitionID to_part) {
    ASSERT(_hg.isBorderNode(hn), "Hypernode " << hn << " is not a border node!");
    DBG(dbg_refinement_kway_fm_move, "moving HN" << hn << " from " << from_part
        << " to " << to_part << " (weight=" << _hg.nodeWeight(hn) << ")");
    _hg.changeNodePart(hn, from_part, to_part);
  }

  PartitionID heaviestPart() const {
    PartitionID heaviest_part = 0;
    for (PartitionID part = 1; part < _config.partition.k; ++part) {
      if (_hg.partWeight(part) > _hg.partWeight(heaviest_part)) {
        heaviest_part = part;
      }
    }
    return heaviest_part;
  }

  void reCalculateHeaviestPartAndItsWeight(PartitionID& heaviest_part,
                                           HypernodeWeight& heaviest_part_weight,
                                           const PartitionID from_part,
                                           const PartitionID to_part) const {
    if (heaviest_part == from_part) {
      heaviest_part = heaviestPart();
      heaviest_part_weight = _hg.partWeight(heaviest_part);
    } else if (_hg.partWeight(to_part) > heaviest_part_weight) {
      heaviest_part = to_part;
      heaviest_part_weight = _hg.partWeight(to_part);
    }
    ASSERT(heaviest_part_weight == _hg.partWeight(heaviestPart()),
           V(heaviest_part) << V(heaviestPart()));
  }

  void reset() {
    _pq.clear();
    _hg.resetHypernodeState();
    _performed_moves.clear();
  }

  Hypergraph& _hg;
  const Configuration& _config;
  KWayRefinementPQ _pq;
  std::vector<RollbackElement> _performed_moves;
  std::vector<HypernodeID> _hns_to_activate;
};
}  // namespace kahypar
