/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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
#include <string>
#include <vector>
#include <set>

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/definitions.h"
#include "kahypar/macros.h"
#include "kahypar/partition/coarsening/rating_tie_breaking_policies.h"
#include "kahypar/partition/preprocessing/louvain.h"
#include "kahypar/partition/preprocessing/quality_measure.h"

namespace kahypar {
class MLCoarsener final : public ICoarsener,
                          private VertexPairCoarsenerBase<>{
 private:
  static const bool dbg_partition_rating = false;

  using Base = VertexPairCoarsenerBase;

  static constexpr HypernodeID kInvalidTarget = std::numeric_limits<HypernodeID>::max();
  static constexpr RatingType kInvalidScore = std::numeric_limits<RatingType>::min();

  struct Rating {
    Rating(HypernodeID trgt, RatingType val, bool is_valid) :
      target(trgt),
      value(val),
      valid(is_valid) { }

    Rating() :
      target(kInvalidTarget),
      value(kInvalidScore),
      valid(false) { }

    Rating(const Rating&) = delete;
    Rating& operator= (const Rating&) = delete;

    Rating(Rating&&) = default;
    Rating& operator= (Rating&&) = delete;

    HypernodeID target;
    RatingType value;
    bool valid;
  };

 public:
  MLCoarsener(Hypergraph& hypergraph, const Configuration& config,
              const HypernodeWeight weight_of_heaviest_node) :
    Base(hypergraph, config, weight_of_heaviest_node),
    _tmp_ratings(_hg.initialNumNodes()), _comm(_hg.initialNumNodes(),0), _louvain(hypergraph,_config) { }

  virtual ~MLCoarsener() { }

  MLCoarsener(const MLCoarsener&) = delete;
  MLCoarsener& operator= (const MLCoarsener&) = delete;

  MLCoarsener(MLCoarsener&&) = delete;
  MLCoarsener& operator= (MLCoarsener&&) = delete;

 private:
  void coarsenImpl(const HypernodeID limit) override final {
    bool ignoreCommunities = false;
    int pass_nr = 0;
    std::vector<HypernodeID> current_hns;
    ds::FastResetFlagArray<> already_matched(_hg.initialNumNodes());
    while (_hg.currentNumNodes() > limit) {
      LOGVAR(pass_nr);
      LOGVAR(_hg.currentNumNodes());
      LOGVAR(_hg.currentNumEdges());
      
      if(_config.preprocessing.use_louvain && !ignoreCommunities) {
        performLouvainCommunityDetection();
      }

      already_matched.reset();
      current_hns.clear();

      const HypernodeID num_hns_before_pass = _hg.currentNumNodes();
      for (const HypernodeID hn : _hg.nodes()) {
        current_hns.push_back(hn);
      }
      Randomize::instance().shuffleVector(current_hns, current_hns.size());
      // std::sort(std::begin(current_hns), std::end(current_hns),
      //           [&](const HypernodeID l, const HypernodeID r) {
      //             return _hg.nodeDegree(l) < _hg.nodeDegree(r);
      //           });
      for (const HypernodeID hn : current_hns) {
        if (_hg.nodeIsEnabled(hn)) {
          const Rating rating = contractionPartner(hn, already_matched);

          if (rating.target != kInvalidTarget) {
            already_matched.set(hn, true);
            already_matched.set(rating.target, true);
            
            performContraction(hn, rating.target);
            if(_config.preprocessing.use_multilevel_louvain) {
                _louvain.contractHypernodes(hn,rating.target);
                size_t N = _hg.initialNumNodes();
                int one_pin_hes_begin = _history.back().one_pin_hes_begin;
                int one_pin_hes_size = _history.back().one_pin_hes_size;
                for(int i = one_pin_hes_begin; i < one_pin_hes_begin+one_pin_hes_size; ++i) {
                    _louvain.contractHypernodes(hn,N+_hypergraph_pruner.removedSingleNodeHyperedges()[i]);
                }
                
                int parallel_hes_begin = _history.back().parallel_hes_begin;
                int parallel_hes_size = _history.back().parallel_hes_size;
                for(int i = parallel_hes_begin; i < parallel_hes_begin+parallel_hes_size; ++i) {
                    _louvain.contractHypernodes(N+_hypergraph_pruner.removedParallelHyperedges()[i].representative_id,
                                                  N+_hypergraph_pruner.removedParallelHyperedges()[i].removed_id);
                }
            }
            
          }

          if (_hg.currentNumNodes() <= limit) {
            break;
          }
        }
      }
      if (num_hns_before_pass == _hg.currentNumNodes()) {
          if(_config.preprocessing.use_louvain && !ignoreCommunities && !_config.preprocessing.only_community_contraction_allowed && _hg.currentNumNodes() >= 2*limit) {
            _comm.assign(_comm.size(),0);
            ignoreCommunities = true;
          }
          else {
            break;
          }
      }

      ++pass_nr;
    }
    //abort();
  }
  
  void performLouvainCommunityDetection() {
      if(_config.preprocessing.use_multilevel_louvain || !_louvain.wasAlreadyExecuted()) {
          HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
          std::set<ClusterID> distinct_comm;
          EdgeWeight quality = _louvain.louvain();   
          for(HypernodeID hn : _hg.nodes()) {
              _comm[hn] = _louvain.clusterID(hn);
            distinct_comm.insert(_comm[hn]);
          }
          HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double> elapsed_seconds = end - start;
          LOG("Louvain-Time: " << elapsed_seconds.count() << "s");
          Stats::instance().addToTotal(_config,"louvainTime",elapsed_seconds.count());
          Stats::instance().addToTotal(_config,"communities",distinct_comm.size()-Stats::instance().get("communities"));
          Stats::instance().addToTotal(_config,"modularity",quality-Stats::instance().get("modularity")); 
      }
  }

  Rating contractionPartner(const HypernodeID u, const ds::FastResetFlagArray<>& already_matched) {
    DBG(dbg_partition_rating, "Calculating rating for HN " << u);
    const HypernodeWeight weight_u = _hg.nodeWeight(u);
    for (const HyperedgeID he : _hg.incidentEdges(u)) {
      ASSERT(_hg.edgeSize(he) > 1, V(he));
      if (_hg.edgeSize(he) <= _config.partition.hyperedge_size_threshold) {
        const RatingType score = static_cast<RatingType>(_hg.edgeWeight(he)) / (_hg.edgeSize(he) - 1);
        for (const HypernodeID v : _hg.pins(he)) {
          if ((v != u && belowThresholdNodeWeight(weight_u, _hg.nodeWeight(v)))) {
            _tmp_ratings[v] += score;
          }
        }
      }
    }
    if (_tmp_ratings.contains(u)) {
      _tmp_ratings.remove(u);
    }

    RatingType max_rating = std::numeric_limits<RatingType>::min();
    HypernodeID target = std::numeric_limits<HypernodeID>::max();
    for (auto it = _tmp_ratings.end() - 1; it >= _tmp_ratings.begin(); --it) {
      const HypernodeID tmp_target = it->key;
      const RatingType tmp_rating = it->value;
      DBG(false, "r(" << u << "," << tmp_target << ")=" << tmp_rating);
      if (_comm[u] == _comm[tmp_target] && acceptRating(tmp_rating, max_rating, target, tmp_target, already_matched)) {
        max_rating = tmp_rating;
        target = tmp_target;
      }
    }

    Rating ret;
    if (max_rating != std::numeric_limits<RatingType>::min()) {
      ASSERT(target != std::numeric_limits<HypernodeID>::max());
      ASSERT(_tmp_ratings[target] == max_rating, V(target));
      ret.value = max_rating;
      ret.target = target;
      ret.valid = true;
    }

    _tmp_ratings.clear();

    ASSERT(!ret.valid || (_hg.partID(u) == _hg.partID(ret.target)));
    DBG(dbg_partition_rating, "rating=(" << ret.value << "," << ret.target << ","
        << ret.valid << ")");
    return ret;
  }

  bool acceptRating(const RatingType tmp, const RatingType max_rating,
                    const HypernodeID old_target, const HypernodeID new_target,
                    const ds::FastResetFlagArray<>& already_matched) const {
    return max_rating < tmp ||
           ((max_rating == tmp) &&
            ((already_matched[old_target] && !already_matched[new_target]) ||
             (already_matched[old_target] && already_matched[new_target] &&
              RandomRatingWins::acceptEqual()) ||
             (!already_matched[old_target] && !already_matched[new_target] &&
              RandomRatingWins::acceptEqual())));
  }

  bool uncoarsenImpl(IRefiner& refiner) override final {
    return doUncoarsen(refiner);
  }

  std::string policyStringImpl() const override final {
    return std::string("");
  }

  bool belowThresholdNodeWeight(const HypernodeWeight weight_u,
                                const HypernodeWeight weight_v) const {
    return weight_v + weight_u <= _config.coarsening.max_allowed_node_weight;
  }

  using Base::_pq;
  using Base::_hg;
  using Base::_config;
  using Base::_history;
  using Base::_hypergraph_pruner;
  ds::SparseMap<HypernodeID, RatingType> _tmp_ratings;
  std::vector<ClusterID> _comm;
  Louvain<Modularity> _louvain;
};
}  // namespace kahypar
