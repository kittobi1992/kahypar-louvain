/***************************************************************************
 *  Copyright (C) 2015 Tobias Heuer <tobias.heuer@gmx.net>
 **************************************************************************/

#ifndef SRC_PARTITION_INITIAL_PARTITIONING_POLICIES_GREEDYQUEUESELECTIONPOLICY_H_
#define SRC_PARTITION_INITIAL_PARTITIONING_POLICIES_GREEDYQUEUESELECTIONPOLICY_H_

#include <limits>

#include "lib/datastructure/FastResetBitVector.h"
#include "lib/datastructure/KWayPriorityQueue.h"
#include "lib/definitions.h"
#include "tools/RandomFunctions.h"

using defs::HypernodeID;
using defs::HyperedgeID;
using defs::PartitionID;
using defs::Hypergraph;
using datastructure::KWayPriorityQueue;


namespace partition {

template<typename GainType = HyperedgeWeight>
struct RoundRobinQueueSelectionPolicy {
  using KWayRefinementPQ = KWayPriorityQueue<HypernodeID, GainType,
                                           std::numeric_limits<GainType>,
                                           ArrayStorage<HypernodeID>, true>;
  //Method returns the part which all hypernodes has to be assigned to before
  //initial partitioning. In experimental results we recognize that it is 
  //desirable to let all hypernodes unassigned before initial partitioning
  //in this method.
  static inline PartitionID getOperatingUnassignedPart() {
    return -1;
  }

  static inline bool nextQueueID(Hypergraph& UNUSED(hg), Configuration& config,
                                 KWayRefinementPQ& _pq, HypernodeID& current_hn, GainType& current_gain,
                                 PartitionID& current_id, bool UNUSED(is_upper_bound_released)) {
    current_id = ((current_id + 1) % config.initial_partitioning.k);
    current_hn = invalid_node;
    current_gain = invalid_gain;
    PartitionID counter = 1;
    while (!_pq.isEnabled(current_id)) {
      if (counter++ == config.initial_partitioning.k) {
        current_id = invalid_part;
        return false;
      }
      current_id = ((current_id + 1) % config.initial_partitioning.k);
    }
    if (current_id != -1) {
      _pq.deleteMaxFromPartition(current_hn, current_gain, current_id);
    }
    return true;
  }

  static const HypernodeID invalid_node = -1;
  static const PartitionID invalid_part = -1;
  static const Gain invalid_gain = std::numeric_limits<Gain>::max();
};

template<typename GainType = HyperedgeWeight>
struct GlobalQueueSelectionPolicy {
    using KWayRefinementPQ = KWayPriorityQueue<HypernodeID, GainType,
                                           std::numeric_limits<GainType>,
                                           ArrayStorage<HypernodeID>, true>;
  //Method returns the part which all hypernodes has to be assigned to before
  //initial partitioning. In experimental results we recognize that it is 
  //desirable to assign all hypernodes to part 1 before initial partitioning
  //in this method.
  static inline PartitionID getOperatingUnassignedPart() {
    return 1;
  }

  static inline bool nextQueueID(Hypergraph& UNUSED(hg), Configuration& config,
                                 KWayRefinementPQ& _pq, HypernodeID& current_hn, GainType& current_gain,
                                 PartitionID& current_id, bool UNUSED(is_upper_bound_released)) {
    current_id = invalid_part;
    current_hn = invalid_node;
    current_gain = invalid_gain;
    bool exist_enabled_pq = false;
    for (PartitionID part = 0; part < config.initial_partitioning.k; ++part) {
      if (_pq.isEnabled(part)) {
        exist_enabled_pq = true;
        break;
      }
    }

    if (exist_enabled_pq) {
      _pq.deleteMax(current_hn, current_gain, current_id);
    }

    ASSERT([&]() {
        if (current_id != -1) {
          for (PartitionID part = 0; part < config.initial_partitioning.k; ++part) {
            if (_pq.isEnabled(part) && _pq.maxKey(part) > current_gain) {
              return false;
            }
          }
        }
        return true;
      } (), "Moving hypernode " << current_hn << " to part "
           << current_id << " isn't the move with maximum gain!");

    return current_id != -1;
  }

  static const HypernodeID invalid_node = -1;
  static const PartitionID invalid_part = -1;
  static const Gain invalid_gain = std::numeric_limits<Gain>::max();
};

template<typename GainType = HyperedgeWeight>
struct SequentialQueueSelectionPolicy {
    using KWayRefinementPQ = KWayPriorityQueue<HypernodeID, GainType,
                                           std::numeric_limits<GainType>,
                                           ArrayStorage<HypernodeID>, true>;
  //Method returns the part which all hypernodes has to be assigned to before
  //initial partitioning. In experimental results we recognize that it is 
  //desirable to assign all hypernodes to part 1 before initial partitioning
  //in this method.
  static inline PartitionID getOperatingUnassignedPart() {
    return 1;
  }

  static inline bool nextQueueID(Hypergraph& hg, Configuration& config,
                                 KWayRefinementPQ& _pq, HypernodeID& current_hn, GainType& current_gain,
                                 PartitionID& current_id, bool is_upper_bound_released) {
    if (!is_upper_bound_released) {
      bool next_part = false;
      if (hg.partWeight(current_id)
          < config.initial_partitioning.upper_allowed_partition_weight[current_id]) {
        _pq.deleteMaxFromPartition(current_hn, current_gain, current_id);

        if (hg.partWeight(current_id) + hg.nodeWeight(current_hn)
            > config.initial_partitioning.upper_allowed_partition_weight[current_id]) {
          _pq.insert(current_hn, current_id, current_gain);
          next_part = true;
        }
      } else {
        next_part = true;
      }

      if (next_part) {
        current_id++;
        if (current_id == config.initial_partitioning.unassigned_part) {
          current_id++;
        }
        if (current_id != config.initial_partitioning.k) {
          ASSERT(_pq.isEnabled(current_id),
                 "PQ " << current_id << " should be enabled!");
          _pq.deleteMaxFromPartition(current_hn, current_gain,
                                     current_id);
        } else {
          current_hn = invalid_node;
          current_gain = invalid_gain;
          current_id = invalid_part;
        }
      }
    } else {
      GlobalQueueSelectionPolicy<GainType>::nextQueueID(hg, config, _pq, current_hn, current_gain,
                                              current_id, is_upper_bound_released);
    }

    return current_id != -1;
  }

  static const HypernodeID invalid_node = -1;
  static const PartitionID invalid_part = -1;
  static const Gain invalid_gain = std::numeric_limits<Gain>::max();
};
}  // namespace partition

#endif  // SRC_PARTITION_INITIAL_PARTITIONING_POLICIES_GREEDYQUEUESELECTIONPOLICY_H_
