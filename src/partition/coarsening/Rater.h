/***************************************************************************
 *  Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_PARTITION_COARSENING_RATER_H_
#define SRC_PARTITION_COARSENING_RATER_H_

#include <limits>
#include <stack>
#include <vector>

#include "lib/datastructure/SparseMap.h"
#include "lib/definitions.h"
#include "lib/macros.h"
#include "partition/Configuration.h"
#include "partition/coarsening/RatingTieBreakingPolicies.h"
#include "lib/datastructure/Neighborhood.h"

using datastructure::SparseMap;
using datastructure::Neighborhood;
using defs::Hypergraph;
using defs::HypernodeID;
using defs::HyperedgeID;
using defs::HypernodeWeight;

namespace partition {
static const bool dbg_partition_rating = false;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
// See Modern C++ Design for the reason why _TiebreakingPolicy has protected non-virtual destructor
template <typename _RatingType, class _TieBreakingPolicy>
class Rater {
 public:
  using RatingType = _RatingType;

 private:
  using TieBreakingPolicy = _TieBreakingPolicy;

  struct HeavyEdgeRating {
    HeavyEdgeRating(HypernodeID trgt, RatingType val, bool is_valid) noexcept :
      target(trgt),
      value(val),
      valid(is_valid) { }

    HeavyEdgeRating() :
      target(std::numeric_limits<HypernodeID>::max()),
      value(std::numeric_limits<RatingType>::min()),
      valid(false) { }

    HeavyEdgeRating(const HeavyEdgeRating&) = delete;
    HeavyEdgeRating& operator= (const HeavyEdgeRating&) = delete;

    HeavyEdgeRating(HeavyEdgeRating&&) = default;
    HeavyEdgeRating& operator= (HeavyEdgeRating&&) = delete;

    HypernodeID target;
    RatingType value;
    bool valid;
  };

 public:
  using Rating = HeavyEdgeRating;

  Rater(Hypergraph& hypergraph, const Configuration& config) noexcept :
    _hg(hypergraph),
    _config(config),
    _tmp_ratings(_hg.initialNumNodes()),
    _neighborhood(hypergraph) { }

  Rater(const Rater&) = delete;
  Rater& operator= (const Rater&) = delete;

  Rater(Rater&&) = delete;
  Rater& operator= (Rater&&) = delete;

  HeavyEdgeRating rate(const HypernodeID u) noexcept {
    DBG(dbg_partition_rating, "Calculating rating for HN " << u);
    const HypernodeWeight weight_u = _hg.nodeWeight(u);
    const PartitionID part_u = _hg.partID(u);
    for (const HyperedgeID he : _hg.incidentEdges(u)) {
      ASSERT(_hg.edgeSize(he) > 1, V(he));
      const RatingType score = static_cast<RatingType>(_hg.edgeWeight(he)) / (_hg.edgeSize(he) - 1);
      for (const HypernodeID v : _hg.pins(he)) {
        if (v != u &&
            belowThresholdNodeWeight(weight_u, _hg.nodeWeight(v)) &&
            (part_u == _hg.partID(v))) {
          _tmp_ratings[v] += score;
        }
      }
    }

    RatingType max_rating = std::numeric_limits<RatingType>::min();
    HypernodeID target = std::numeric_limits<HypernodeID>::max();
    for (auto it = _tmp_ratings.crbegin(); it != _tmp_ratings.crend(); ++it) {
      const HypernodeID tmp_target = it->key;
      const RatingType tmp = it->value /
                             (weight_u * _hg.nodeWeight(tmp_target));
      DBG(false, "r(" << u << "," << tmp_target << ")=" << tmp);
      double jaccard_index = 0.0;
      double tmp_index = _neighborhood.jaccardIndex(u, tmp_target);
      if (acceptRating(tmp, max_rating) && jaccard_index < tmp_index) {
        max_rating = tmp;
        target = tmp_target;
        jaccard_index = tmp_index;
      }
    }
    _tmp_ratings.clear();
    HeavyEdgeRating ret;
    if (max_rating != std::numeric_limits<RatingType>::min()) {
      ASSERT(target != std::numeric_limits<HypernodeID>::max(), "invalid contraction target");
      ret.value = max_rating;
      ret.target = target;
      ret.valid = true;
    }
    ASSERT([&]() {
        bool flag = true;
        if (ret.valid && (_hg.partID(u) != _hg.partID(ret.target))) {
          flag = false;
        }
        return flag;
      } (), "Representative " << u << " & contraction target " << ret.target
           << " are in different parts!");
    DBG(dbg_partition_rating, "rating=(" << ret.value << "," << ret.target << ","
        << ret.valid << ")");
    return ret;
  }

  HypernodeWeight thresholdNodeWeight() const noexcept {
    return _config.coarsening.max_allowed_node_weight;
  }

  Neighborhood& neighborhood() {
    return _neighborhood;
  }

 private:
  bool belowThresholdNodeWeight(const HypernodeWeight weight_u,
                                const HypernodeWeight weight_v) const noexcept {
    return weight_v + weight_u <= _config.coarsening.max_allowed_node_weight;
  }

  bool acceptRating(const RatingType tmp, const RatingType max_rating) const noexcept {
    return max_rating < tmp || (max_rating == tmp);
  }

  Hypergraph& _hg;
  const Configuration& _config;
  SparseMap<HypernodeID, RatingType> _tmp_ratings;
  std::vector<RatingType> _tmp_ratings;
  Neighborhood _neighborhood;
};
#pragma GCC diagnostic pop
}  // namespace partition

#endif  // SRC_PARTITION_COARSENING_RATER_H_
