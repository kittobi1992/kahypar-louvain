/***************************************************************************
 *  Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_PARTITION_REFINEMENT_FMREFINERBASE_H_
#define SRC_PARTITION_REFINEMENT_FMREFINERBASE_H_

#include <limits>

#include "lib/definitions.h"
#include "partition/Configuration.h"
#include "partition/refinement/DeltaGainCollector.h"

using defs::Hypergraph;
using defs::HypernodeID;
using defs::HyperedgeID;

using partition::DeltaGainCollector;

namespace partition {
static const bool dbg_refinement_fm_border_node_check = false;
static const bool dbg_refinement_kway_fm_move = false;

class FMRefinerBase {
 public:
  using Gain = HyperedgeWeight;

 protected:
  static constexpr HypernodeID kInvalidHN = std::numeric_limits<HypernodeID>::max();
  static constexpr Gain kInvalidGain = std::numeric_limits<Gain>::min();
  static constexpr HyperedgeWeight kInvalidDecrease = std::numeric_limits<PartitionID>::min();

  FMRefinerBase(Hypergraph& hypergraph, const Configuration& config) noexcept :
    _hg(hypergraph),
    _config(config),
    _hypernode_connected_to_part(_hg.initialNumNodes(),_config.partition.k) { }

  ~FMRefinerBase() { }

  FMRefinerBase(const FMRefinerBase&) = delete;
  FMRefinerBase& operator= (const FMRefinerBase&) = delete;

  FMRefinerBase(FMRefinerBase&&) = delete;
  FMRefinerBase& operator= (FMRefinerBase&&) = delete;

  bool hypernodeIsConnectedToPart(const HypernodeID pin, const PartitionID part) noexcept {
    std::pair<HypernodeID,PartitionID> key = std::make_pair(pin,part);
    if(_hypernode_connected_to_part.isValidEntry(key)) {
      return _hypernode_connected_to_part.get(key);
    }
    for (const HyperedgeID he : _hg.incidentEdges(pin)) {
      if (_hg.pinCountInPart(he, part) > 0) {
	_hypernode_connected_to_part.add(key,true);
        return true;
      }
    }
    _hypernode_connected_to_part.add(key,false);
    return false;
  }

  bool moveIsFeasible(const HypernodeID max_gain_node, const PartitionID from_part,
                      const PartitionID to_part) const noexcept {
    ASSERT(_config.partition.mode == Mode::direct_kway,
           "Method should only be called in direct partitioning");
    return (_hg.partWeight(to_part) + _hg.nodeWeight(max_gain_node)
            <= _config.partition.max_part_weights[0]) && (_hg.partSize(from_part) - 1 != 0);
  }

  void moveHypernode(const HypernodeID hn, const PartitionID from_part,
                     const PartitionID to_part) noexcept {
    ASSERT(_hg.isBorderNode(hn), "Hypernode " << hn << " is not a border node!");
    DBG(dbg_refinement_kway_fm_move, "moving HN" << hn << " from " << from_part
        << " to " << to_part << " (weight=" << _hg.nodeWeight(hn) << ")");
    _hypernode_connected_to_part.clear();
    _hg.changeNodePart(hn, from_part, to_part);
  }

  PartitionID heaviestPart() const noexcept {
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
                                           const PartitionID to_part) const noexcept {
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

  Hypergraph& _hg;
  const Configuration& _config;
  DeltaGainCollector<bool> _hypernode_connected_to_part;
};
}  // namespace partition
#endif  // SRC_PARTITION_REFINEMENT_FMREFINERBASE_H_
